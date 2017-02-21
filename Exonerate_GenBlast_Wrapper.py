#!/usr/bin/env python3
import re
import os
from collections import defaultdict, namedtuple
import tempfile as tmp
from shared_code_box import run_cmd, write_to_tempfile, hash_fasta, get_phmm_score, tempdir


########################################################################################################################
# Prediction Object and related functions
########################################################################################################################

def isolate_overlapping_predictions(pred_obj_list):
    passed = []
    previous_interactions = set()
    for pred_obj in pred_obj_list:
        region_owner = pred_obj
        if pred_obj not in previous_interactions:
            for next_pred_obj in pred_obj_list:
                if next_pred_obj not in previous_interactions and pred_obj.check_overlap(next_pred_obj):
                    previous_interactions.update([pred_obj, next_pred_obj])     # can still be same cluster but overlapping on different strands!
                    competitor = [next_pred_obj, region_owner]
                    region_owner = sorted(competitor, key=lambda x: (x.score, x.aln_score), reverse=True)[0]
            passed.append(region_owner)
    return previous_interactions.difference(set(passed)), passed


class PredictionObject:
    def __init__(self, region, score, cluster, cutoff, length_range):
        self.cluster = cluster
        self.score = score
        self.cutoff = cutoff
        self.region = region
        self.length_range = length_range
        self.strand = self.region.strand
        self.contig = self.region.contig
        self.DNA = None
        self.protein = None
        self.gff = None
        self.blast_location = None
        self.gene_start = None
        self.gene_end = None
        self.gene_length = None
        self.aln_score = None
        self.fragmented = None
        self.geneID = None

    def infer_data_from_exonerate_obj(self, exonerate_obj, key_tuple):
        if type(exonerate_obj) == ExonerateObject:
            self.aln_score = key_tuple.score
            self.blast_location = (int(self.region.s_start), int(self.region.s_end))
            self.DNA = exonerate_obj.target_dna[key_tuple]
            self.protein = exonerate_obj.target_prot[key_tuple]
            self.gff = write_exonerate_gff(exonerate_obj.gff[key_tuple], self.blast_location, self.strand)
            self.strand = self.gff[0].split("\t")[6]
            self.gene_start, self.gene_end = [int(x) for x in self.gff[0].split("\t")[3:5]]
            self.gene_length = abs(self.gene_end - self.gene_start)
            self.geneID = self.gff[0].split("\t")[-1].split('''"''')[1]
        else:
            return None

    def fragmentation_check(self):
        if float(self.length_range[0]) <= len(self.protein) <= float(self.length_range[1]):
            self.fragmented = False
        else:
            self.fragmented = True
        return self.fragmented

    def check_overlap(self, other):
        if self != other and self.contig == other.contig:
            coordinates = [other.gene_start, other.gene_end, self.gene_start, self.gene_end]
            if (other.gene_length + self.gene_length) >= (max(coordinates) - min(coordinates)):
                if self.strand == other.strand:
                    return True
                else:
                    return False
        return False


########################################################################################################################
# Exonerate Wrapper
########################################################################################################################


aa3_to_1_coding_dict = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
     'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N',
     'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W', '***': 'X',
     'Ala': 'A', 'Val':'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M','Xaa': 'X', 'Unk': 'X'}


def aacode_3to1(seq):
    """Turn a three letter protein into a one letter protein.
The 3 letter code can be upper, lower, or any mix of cases
The seq input length should be a factor of 3 or else results
in an error
"""
    if len(seq) % 3 == 0:
        single_seq = []
        for i in range(0, len(seq), 3):
            single_seq.append(aa3_to_1_coding_dict.get(seq[i:i+3]))
        return "".join(single_seq)
    else:
        raise False


def kill_char(string, n):
    if n > len(string):
        raise IndexError
    begin = string[:n]
    end = string[n+1:]
    return begin + end


def find_hashes(protein_string):
    for idx, sign in enumerate(protein_string):
        if sign == "#":
            yield idx


def clear_hashed_bases(protein_string, dna_string):
    """removes nucleotides which are opposite to #'s in the protein sequence (non-coding-nucleotides)"""
    if "#" in protein_string:
        for idx in sorted(find_hashes(protein_string), reverse=True):
            dna_string = kill_char(dna_string, idx)
    return dna_string


def remove_non_letter_signs(string):
    """removes all non-alphabetic letters from string"""
    regex = re.compile('[^a-zA-Z*]')
    return regex.sub("", string)


def remove_lower(text_string):
    return re.sub('[a-z]', '', text_string)


def all_proteins_to_fasta_string(exonerate_obj):
    """returns a fasta string containing all predicted protein sequences, with query name as header."""
    fasta_list = []
    for hit_tuple in exonerate_obj.query_prot:
        fasta_list.append(">" + hit_tuple.query + ";{}".format(hit_tuple.idx))
        fasta_list.append(exonerate_obj.target_prot[hit_tuple])
    return "\n".join(fasta_list)


def correct_gene_position(strand, blast_location, exonerate_location):
    if strand == "-":
        gene_start = max(blast_location) - (min(exonerate_location) - 1)
        gene_end = max(blast_location) - (max(exonerate_location) - 1)
    elif strand == "+":
        gene_start = min(blast_location) + (min(exonerate_location) - 1)
        gene_end = min(blast_location) + (max(exonerate_location) - 1)
    else:
        print("[!]\t ERROR: {} can not be identified as + or - strand".format(strand))
        return None
    return min([str(gene_start), str(gene_end)]), max([str(gene_start), str(gene_end)])


def write_exonerate_gff(gff_list, off_set_tuple, strand):
    last_phase = 0
    new_gff_list = []
    cds_count = 1
    gene_name, rna_name = "", ""
    for line in gff_list:
        line = line.split("\t")
        line[3], line[4] = correct_gene_position(line[6], off_set_tuple, (int(line[3]), int(line[4])))
        if "cds" in line[2]:
            current_phase = (abs(int(line[4]) - int(line[3])) + 1 + last_phase) % 3
            line[7] = str(last_phase)
            if current_phase >= 3:
                last_phase = 0
            else:
                last_phase = current_phase
            line[-1] = "{}; {}; exon_number '{}';".format(gene_name, rna_name, str(cds_count))
            cds_count += 1
        elif "gene" in line[2]:
            name = line[-1].split(";")[1].split(" ")[2].strip() + "_{}".format(str(line[3]))
            gene_name = "gene_id " + '''"{}"'''.format(name)
            rna_name = "transcript_id " + '''"{}"'''.format(name + ".t")
            line[2] = "transcript"
            line[-1] = "{}; {};".format(gene_name, rna_name)
        else:
            continue
        new_gff_list.append("\t".join(line))
    return new_gff_list


def make_exonerate_command(model, query_file, region_file):
    cmd = "exonerate {} --softmaskquery no -Q protein -T dna  " \
      "--softmasktarget no --showvulgar no --minintron 20 --maxintron 50000 " \
      "--showalignment yes --showtargetgff yes -q {} -t {}".format(model, query_file, region_file) #--minintron 20 --maxintron 50000
    return cmd


def get_exonerate_object(output_path, command):
    line_count = 0
    with open(output_path, "w") as ex:
        for line in run_cmd(command=command, wait=False):
            ex.write(line)
            line_count += 1
    if line_count < 10:
        return None
    else:
        return ExonerateObject(output_path)


def run_exonerate(mode_string, name, directory, region, query):    # sometimes exhaustive results in core dump
    cmd = make_exonerate_command(mode_string, query, region)     # "-m p2g:b -E yes" is the best exhaustive command
    out_file = os.path.join(directory, name)
    exonerate_obj = get_exonerate_object(out_file, cmd)
    if exonerate_obj is None and "-E yes" in mode_string:
        print("coredump")
        cmd = make_exonerate_command("-m p2g -E no", query, region)
        exonerate_obj = get_exonerate_object(out_file, cmd)
    return exonerate_obj


class ExonerateObject:

    def __init__(self, exonerate_file):
        self.path = os.path.realpath(exonerate_file)
        self.query_prot = defaultdict(list)
        self.target_dna = defaultdict(list)
        self.target_prot = defaultdict(list)
        self.gff = defaultdict(list)
        self.exonerate_processor(exonerate_file)

    def exonerate_processor(self, exonerate_file):
        ''' protein against 1 or more genome regions (targets),
        file with aln and/or gff output '''
        del_intron = lambda seq: seq.replace("TargetIntron", "")                  # no lambda
        prep_line = lambda line: line.strip("\n").split(": ")[1].strip(" ")
        next_line = lambda: prep_line(next(ex))
        next_block = lambda size: [next_line() for x in range(0, size)]

        read_flag = 0
        Region = namedtuple('Region', 'query, qrange, trange, target, score, idx')
        with open(exonerate_file, "r") as ex:
            count = 0
            for line in ex:
                if line.strip():
                    if read_flag == 0:
                        if line.startswith("C4 Alignment"):
                            read_flag = 1
                    elif read_flag == 1:
                        count += 1
                        query, target, model, score = next_block(4)
                        qrange, trange = [tuple(x.split(" -> ")) for x in next_block(2)]
                        header_tuple = Region(query=query, qrange=qrange, trange=trange, target=target, score=score, idx=count)
                        read_flag = 2
                    elif read_flag == 2:
                        if not line.startswith("#"):
                            self.query_prot[header_tuple].append(del_intron(remove_non_letter_signs(line)))
                            next(ex)
                            target_prot = next(ex)
                            target_dna = clear_hashed_bases(target_prot, next(ex))
                            self.target_prot[header_tuple].append(remove_non_letter_signs(target_prot))
                            self.target_dna[header_tuple].append(remove_lower(remove_non_letter_signs(target_dna)))
                        else:
                            self.query_prot[header_tuple] = "".join(self.query_prot[header_tuple])
                            self.target_dna[header_tuple] = "".join(self.target_dna[header_tuple])
                            self.target_prot[header_tuple] = aacode_3to1("".join(self.target_prot[header_tuple]))
                            if "GFF" in line:
                                read_flag = 3
                            else:
                                read_flag = 0
                    elif read_flag == 3:
                        if not line.startswith("#"):
                            line = line.strip("\n")
                            self.gff[header_tuple].append(line)
                        elif "END OF GFF DUMP" in line:
                            read_flag = 0
                    else:
                        pass
            ex.seek(0)


def markov_model_scoring2(fasta_string, hmm):
    if hmm:
        with tmp.NamedTemporaryFile() as ex_file:
            write_to_tempfile(ex_file.name, fasta_string)
            score_hash = get_phmm_score(hmm, ex_file.name)
            if score_hash:
                return score_hash
    return None


def best_exonerate_prediction2(region_fasta, query_fasta, dir_path, hmm):
    with tmp.NamedTemporaryFile() as reg_file:
        write_to_tempfile(reg_file.name, region_fasta)
        ex_obj = run_exonerate("-m p2g -E no", "{}.exon_p2g", dir_path, region_fasta, query_fasta)
        if ex_obj:
            all_proteins = all_proteins_to_fasta_string(ex_obj)
            TP_scores = markov_model_scoring2(all_proteins, hmm)
            if TP_scores:
                print(TP_scores)
                max_score = max(list(TP_scores.values()))
                max_score_header = set([header.split(";")[0] for header, score in TP_scores.items() if score >= max_score*0.90])
                fasta_hash = hash_fasta(query_fasta) # do that before
                with tmp.NamedTemporaryFile() as q_file:
                    max_val_fasta = "\n".join(["{}\n{}".format(header, fasta_hash[header]) for header in max_score_header])
                    write_to_tempfile(q_file.name, max_val_fasta)
                    ex_name = "{}.exon".format(q_file.name)
                    ex_obj = run_exonerate("-m p2g -E no ", ex_name, dir_path, region_fasta, q_file.name)
                    return ex_obj
    return None


if __name__ == "__main__":
    with tempdir() as tmp_dir:
        f22_region = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/eef_data/F226Dparalog_T15D6.2_region.fa"
        b0205_region = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/eef_data/c34_region.fa"
        #b0205_protein = "/home/jgravemeyer/Dropbox/MSc_project/data/accuracy_ortho_group_files_5.0/C34F11.1_Inf5.0_OG0004430.fa"
        b0205_protein = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/eef_data/C34F11.1_Inf5.0_OG0004430.fa"
        f22_protein_fasta = "/home/jgravemeyer/Dropbox/MSc_project/data/accuracy_ortho_group_files_5.0/F22D6.12_Inf5.0_OG0000412.fa"
        hmm_f22 = "/home/jgravemeyer/Dropbox/MSc_project/data/out_accuracy_inf5/F22D6.12_Inf5.0_OG0000412.hmm"
        hmm_b0205 = "/home/jgravemeyer/Dropbox/MSc_project/data/out_accuracy_inf5/C34F11.1_Inf5.0_OG0004430.hmm"
    #    test = run_exonerate("-m p2g -E yes", "F22test_exonerate.out", "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/eef_data",
    #                  "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/eef_data/F226Dparalog_T15D6.2_region.fa",
    #                  "/home/jgravemeyer/Dropbox/MSc_project/data/accuracy_ortho_group_files_5.0/F22D6.12_Inf5.0_OG0000412.fa")
        print("########################################################################################################################")

    test2 = run_exonerate("-m p2g -E yes", "F22test_exonerate.out", "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/eef_data",
                  "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/eef_data/F226Dparalog_T15D6.2_region.fa",
                  "/home/jgravemeyer/programs/scipio-1.4/T15D6_F22paralog.fa")


    for predtuple in test2.target_dna:
        print("\n".join(write_exonerate_gff(test2.gff[predtuple], (100000,15000), "+")))

        '''
        ex_obj_f22 = best_exonerate_prediction2(f22_region, f22_protein_fasta, tmp_dir, hmm_f22)
        region_offset = {"all":[12373740, 12384000], "first": [12376000, 12378905], "second": [12378905, 12381250], "b0205": [10704500, 10715000]}
        all_proteins = all_proteins_to_fasta_string(ex_obj_f22)
        TP_scores = markov_model_scoring2(all_proteins, hmm_f22)
        Region = namedtuple('Region', 'query, s_start, s_end, strand, contig')
        region_f22 = Region(query='f22', s_start=12373740, s_end=12384000, strand="-", contig="I")
        f22_pred_obj_list = []
        for pred_tuple in ex_obj_f22.target_prot:
            score_id = ">" + pred_tuple.query + ";{}".format(pred_tuple.idx)
            pred_obj = PredictionObject(region_f22, TP_scores[score_id], pred_tuple.query, "cutoff", "lenghtrange")
            pred_obj.infer_data_from_exonerate_obj(ex_obj_f22, pred_tuple)
            f22_pred_obj_list.append(pred_obj)
        for single_obj in sorted(f22_pred_obj_list, key=lambda x: (x.score, x.aln_score), reverse=True):
            print(single_obj.gene_start, single_obj.gene_end, single_obj.cluster, single_obj.score, len(single_obj.protein))
        print("finally")
        filtered, passed = isolate_overlapping_predictions(sorted(f22_pred_obj_list, key=lambda x: x.score, reverse=True))
        for valid in passed:
            print(valid.gff)


        for headertuple in ex_obj_f22.target_dna:
            print(headertuple)
            offset = correct_gene_position("-", region_offset["all"], [int(headertuple.trange[0]), int(headertuple.trange[1])])
            offset2 = [int(x) for x in sorted(offset)]
            print(len(ex_obj_f22.target_prot[headertuple]), TP_scores[">" + headertuple.query + ";{}".format(headertuple.idx)])
            offset_range = set(range(offset2[0], offset2[1]))
            first_gene = (len(offset_range.intersection(range(12376904, 12378674)))) / 1770
            second_gene = (len(offset_range.intersection(range(12379335, 12381032)))) / 1697
            print(first_gene, "first gene")
            print(second_gene, "second gene")
        print("#"*30)


        print("#"*30)
        ex_obj_b20 = best_exonerate_prediction2(b0205_region, b0205_protein, tmp_dir, hmm_b0205)
        region_offset = {"all": [12373740, 12384000], "first": [12376000, 12378905], "second": [12378905, 12381250], "b0205": [10704500, 10715000], "c34": [9948644,9964239], "zk353": [8405000, 8415000]}
        all_proteins = all_proteins_to_fasta_string(ex_obj_b20)
        TP_scores = markov_model_scoring2(all_proteins, hmm_b0205)
        Region = namedtuple('Region', 'query, s_start, s_end, strand, contig')
        region_b20 = Region(query='c34', s_start=9948644, s_end=9964239, strand="-", contig="I")
        b20_pred_obj_list = []
        for pred_tuple in ex_obj_b20.target_prot:
            score_id = ">" + pred_tuple.query + ";{}".format(pred_tuple.idx)
            if score_id in TP_scores:
                pred_obj = PredictionObject(region_b20, TP_scores[score_id], pred_tuple.query, "cutoff", "lenghtrange")
                pred_obj.infer_data_from_exonerate_obj(ex_obj_b20, pred_tuple)
                b20_pred_obj_list.append(pred_obj)
        filtered, passed = isolate_overlapping_predictions(b20_pred_obj_list)
        for single_obj in sorted(b20_pred_obj_list, key=lambda x: (x.score, x.aln_score), reverse=True):
            print(single_obj.gene_start, single_obj.gene_end, single_obj.cluster, single_obj.score, len(single_obj.protein), single_obj.aln_score)
            #print("\n".join(single_obj.gff))
        print("finally")
        for valid in passed:
            print(valid.gene_start, valid.gene_end, valid.cluster, valid.score, len(valid.protein), valid.strand)
            print("\n".join(valid.gff))
        for removed in filtered:
            print(removed.gene_start, removed.gene_end)


        for headertuple in ex_obj_b20.target_dna:
            print(headertuple)
            offset = correct_gene_position("-", region_offset["W02D9"], [int(headertuple.trange[0]), int(headertuple.trange[1])])
            offset2 = [int(x) for x in sorted(offset)]
            print(len(ex_obj_b20.target_prot[headertuple]), "lenght", TP_scores[">" + headertuple.query + ";{}".format(headertuple.idx)])
            offset_range = set(range(offset2[0], offset2[1]))
            first_gene = (len(offset_range.intersection(range(12557000,12564481)))) / 451
            print(first_gene, "first gene")
        '''

'''
########################################################################################################################
9953552 9955716 CTROP.Csp11.Scaffold629.g10913.t1 109 401 -
I	exonerate:protein2genome:local	gene	9953552	9955716	546	-	.	gene_id 1 ; sequence CTROP.Csp11.Scaffold629.g10913.t1 ; gene_orientation +
I	exonerate:protein2genome:local	mRNA	9953552	9955716	546	-	.	gene_id 1 ; sequence CTROP.Csp11.Scaffold629.g10913.t1 ; gene_orientation +
I	exonerate:protein2genome:local	cds	9953552	9953948	.	-	0
I	exonerate:protein2genome:local	exon	9953552	9953948	.	-	.	insertions 3 ; deletions 0
I	exonerate:protein2genome:local	intron	9953949	9954336	.	-	.	intron_id 1
I	exonerate:protein2genome:local	cds	9954337	9954403	.	-	1
I	exonerate:protein2genome:local	exon	9954337	9954403	.	-	.	insertions 0 ; deletions 0
I	exonerate:protein2genome:local	intron	9954404	9954460	.	-	.	intron_id 2
I	exonerate:protein2genome:local	cds	9954461	9954534	.	-	2
I	exonerate:protein2genome:local	exon	9954461	9954534	.	-	.	insertions 0 ; deletions 0
I	exonerate:protein2genome:local	intron	9954535	9954564	.	-	.	intron_id 3
I	exonerate:protein2genome:local	cds	9954565	9954618	.	-	1
I	exonerate:protein2genome:local	exon	9954565	9954618	.	-	.	insertions 0 ; deletions 0
I	exonerate:protein2genome:local	intron	9954619	9954696	.	-	.	intron_id 4
I	exonerate:protein2genome:local	cds	9954697	9954900	.	-	1
I	exonerate:protein2genome:local	exon	9954697	9954900	.	-	.	insertions 0 ; deletions 0
I	exonerate:protein2genome:local	intron	9954901	9954927	.	-	.	intron_id 5
I	exonerate:protein2genome:local	cds	9954928	9955090	.	-	1
I	exonerate:protein2genome:local	exon	9954928	9955090	.	-	.	insertions 0 ; deletions 1
I	exonerate:protein2genome:local	intron	9955091	9955223	.	-	.	intron_id 6
I	exonerate:protein2genome:local	cds	9955224	9955316	.	-	2
I	exonerate:protein2genome:local	exon	9955224	9955316	.	-	.	insertions 0 ; deletions 0
I	exonerate:protein2genome:local	intron	9955317	9955565	.	-	.	intron_id 7
I	exonerate:protein2genome:local	cds	9955566	9955716	.	-	2
I	exonerate:protein2genome:local	exon	9955566	9955716	.	-	.	insertions 0 ; deletions 1

Process finished with exit code 0
'''
