#!/usr/bin/env python3

'''
Usage: use_models.py                          -i <DIR> -g <FILE> [-c <INT>] [-o <DIR>] [--keep] [--verbose] [--frag] [--quick]

    Options:
        -h, --help                            show this screen.

        General
        -i, --use_models_input <DIR>          Directory with results from build_models.py
        -g, --genome <FILE|LIST>              Single genome fasta-file or a many lines ".txt"-file in style of: name=path/to/genome
        Optional
        -c, --coverage_filer <INT>            Minimal aligned length of a Blast query to the target genome in % (used to filter Blast hits)
        -o, --out_dir <DIR>                   Directory for the output files and folders (default: same as input directory)
        --quick                               Omits exhaustive prediction refinements through and additional exonerate -E yes run (faster but less accurate)
        --frag                                If enabled, a length filter will be applied to remove potentially fragmented predictions
        --keep                                Keeps intermediate files (Blast output, merged regions, exonerate output)
        --verbose                             Prints progress details to the screen

'''

import os
import sys
import logging
import tempfile as tmp
from docopt import docopt
from shared_code_box import tempdir, check_programs, get_phmm_score, write_to_tempfile, get_outdir, hash_fasta
from Exonerate_GenBlast_Wrapper import run_exonerate, all_proteins_to_fasta_string,\
    isolate_overlapping_predictions, PredictionObject
from Blast_wrapper import run_tblastn, make_blast_db
from collections import defaultdict

########################################################################################################################
# Global Functions
########################################################################################################################
coverage_min = 30
data_base = None
gene_ps_results = None
out_dir = None
keep = None
verbose = None
genome = None
frag = None
quick = None
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logger_blast_region = logging.getLogger("BLAST")
logger_prediction = logging.getLogger("Exonerate")
logger_validate = logging.getLogger("Filtering")


def format_genome_hash(argument):
    genome_hash = {}
    error_list = []
    if args['--genome'].split(".")[-1] == "txt":
        with open(args['--genome']) as g_file:
            for genome_line in g_file:
                try:
                    genome_line = genome_line.strip("\n").split("=")
                    if os.path.exists(genome_line[1]):
                        genome_hash[genome_line[1]] = genome_line[0]
                    else:
                        error_list.append("[!]\t ERROR: {} does not exist or was specified wrongly".format(genome_line[0]))
                except IndexError:
                    error_list.append("[!]\t ERROR: genome file is not a TAB-SEPARATED 2-column file".format(genome_line))
    else:
        try:
            genome_hash[os.path.abspath(args["--genome"])] = os.path.split(args["--genome"])[-1].split(".")[0]
        except IndexError:
            genome_hash[os.path.abspath(args["--genome"])] = os.path.split(args["--genome"])
    return error_list, genome_hash


def check_arguments(args):
    global coverage_min, out_dir, gene_ps_results, keep, verbose, genome, frag, quick
    gene_ps_results = os.path.abspath(args['--use_models_input'])
    keep = args['--keep']
    verbose = args['--verbose']
    frag = args['--frag']
    quick = args['--quick']
    error_list = []
    if not os.path.exists(gene_ps_results):
        error_list.append("[!]\t ERROR: input directory: {} does not exist".format(gene_ps_results))
    if not os.path.isdir(gene_ps_results):
        error_list.append("[!]\t ERROR: please specify a DIRECTORY as input, {} is not a directory".format(gene_ps_results))
    if args['--coverage_filer']:
        try:
            coverage_min = int(args['--coverage_filer'])
        except ValueError:
            error_list.append("[!]\t ERROR: coverage_min needs integer; '{}' is not an integer".format(coverage_min))
    if args["--out_dir"]:
        if os.path.isdir(args["--out_dir"]):
            out_dir = get_outdir(args["--out_dir"])
        else:
            out_dir = get_outdir("/".join(gene_ps_results.split("/")[:-1]), add_dir=args["--out_dir"])
    else:
        out_dir = get_outdir("/".join(gene_ps_results.split("/")[:-1]), add_dir="Predictions")
    errors, genome_hash = format_genome_hash(args)
    error_list.extend(errors)
    if error_list:
        print("[!] {} - ARGUMENT ERRORS\n".format(len(error_list)))
        print("\n".join(error_list))
        sys.exit()
    return genome_hash


def coverage_filter(area):
    if area.chunk_cov == 100 and area.query_cov < 20:
        return False
    elif area.chunk_cov > coverage_min:
        return True


def convert_summary_to_string(summary_array):
    summary_list = []
    row_format = "\t{}" * (len(summary_array[0]) + 1)     # length of header
    for entry in summary_array:
        summary_list.append(row_format.format("", *entry))
    return "\n".join(summary_list)


def markov_model_scoring(fasta_string, hmm):
    if hmm:
        with tmp.NamedTemporaryFile() as ex_file:
            write_to_tempfile(ex_file.name, fasta_string)
            score_hash = get_phmm_score(hmm, ex_file.name)
            if score_hash:
                return score_hash
    return None


def write_merged_region_to_intermediate(blast_ob):
    results_list = [
        "# Merging Distance: {}, Flanking Distance {}".format(blast_ob.merging_distance, blast_ob.flanking_distance),
        "# Fields: contig, subject, subject_start, subject_end, strand, chunk_coverage, total_coverage, query_length"]
    for contig in blast_ob.inferred_regions:
        for query in blast_ob.inferred_regions[contig]:
            for region in blast_ob.inferred_regions[contig][query]:
                results_list.append("\t".join([region.contig, query, str(region.s_start), str(region.s_end), str(region.strand), str(region.chunk_cov), str(region.query_cov), str(region.q_len)]))
    return "\n".join(results_list)


def find_best_exonerate_result(region_tuple, region_fasta, group, cluster, dir_path):
    '''first aligns all proteins in p2g mode, scores them against the hmm and extracts the highest scoring predictions.
    Then the corresponding actual proteins of those are re-aligned in p2g -E yes mode. Returned exonerate object,
    holds thoses (few-best) predictions -> needs for further filtering.'''
    query_fasta = data_base.group_by_cluster_to_fasta_file[group][cluster]
    query_hash = data_base.group_by_cluster_to_fasta_hash[group][cluster]
    hmm = data_base.group_by_cluster_to_hmm[group][cluster]
    cutoff = data_base.group_by_cluster_to_score_cutoff[group][cluster]
    length_range = data_base.group_by_cluster_to_length_range[group][cluster]
    with tmp.NamedTemporaryFile() as reg_file:
        write_to_tempfile(reg_file.name, region_fasta)
        ex_obj = run_exonerate("-m p2g -E no", "{}.exon_p2g", dir_path, reg_file.name, query_fasta)
        try:
            TP_scores = markov_model_scoring(all_proteins_to_fasta_string(ex_obj), hmm)
            if not quick:
                max_score = max(list(TP_scores.values()))
                max_header = set([header.split(";")[0] for header, score in TP_scores.items() if score == max_score])
                with tmp.NamedTemporaryFile() as q_file:
                    max_val_fasta = "\n".join(["{}\n{}".format(header, query_hash[header]) for header in max_header])
                    write_to_tempfile(q_file.name, max_val_fasta)
                    ex_obj = run_exonerate("-m p2g -E yes", "{}.exon".format(q_file.name), dir_path, reg_file.name, q_file.name)
                    TP_scores = markov_model_scoring(all_proteins_to_fasta_string(ex_obj), hmm)
            best_pred = []
            for key_tuple in ex_obj.target_dna:
                score_id = ">" + key_tuple.query + ";{}".format(key_tuple.idx)
                pred_obj = PredictionObject(region_tuple, TP_scores[score_id], cluster, cutoff, length_range)
                pred_obj.infer_data_from_exonerate_obj(ex_obj, key_tuple)
                best_pred.append(pred_obj)
            filtered, passed = isolate_overlapping_predictions(sorted(best_pred, key=lambda x: x.score, reverse=True))
            return ex_obj, passed
        except (AttributeError, TypeError, KeyError):
            return None, None


def prediction_filter(pred_obj, TN_HMM):
    status = False
    TP_score = pred_obj.score
    if TP_score and TP_score >= pred_obj.cutoff:
        if TN_HMM:
            protein_fasta = ">{}\n{}".format(pred_obj.cluster, pred_obj.protein)
            try:
                TN_score = list(markov_model_scoring(protein_fasta, TN_HMM).values())[0]
                if TN_score < TP_score:
                    status = True
            except AttributeError:
                pass
        else:
            status = True
    if frag and status is True:
        if pred_obj.fragmentation_check() is True:
            status = "FRAG"
    return status


########################################################################################################################
# Class: DataProvider (reads in and processes Gene PS files
########################################################################################################################

class DataProviderObject:

    def __init__(self, input_data):
        self.gene_ps_results = input_data
        self.group_names = []
        self.group_to_original_path = {}
        self.group_to_group_size = {}
        self.group_to_consensus_file = {}
        self.group_by_cluster_to_hmm = defaultdict(dict)
        self.group_by_cluster_to_TN_hmm = defaultdict(dict)
        self.group_by_cluster_to_fasta_file = defaultdict(dict)
        self.group_by_cluster_to_fasta_hash = defaultdict(dict)
        self.group_by_cluster_to_length_range = defaultdict(dict)
        self.group_by_cluster_to_score_cutoff = defaultdict(dict)
        self.cluster_scope = self.load_data_and_initialize_global_variables()

    def check_loaded_data(self, cluster_count, input_scope, error_list):
        if cluster_count != input_scope:
            error_list.append("\t[!] Supposed number of cluster {} != counted amount of cluster {}\n".format(cluster_count, input_scope))
        if cluster_count == 0:
            error_list.append("\t[!] No GenePS files found - please check your input directory in terms of path and file style\n")
        if error_list:
            print("\n".join(error_list))
            sys.exit()
        else:
            return input_scope

    def validate_path_files(self, cluster):
        error_list = []
        for ending in [".hmm", ".fasta"]:
            if os.path.exists(os.path.join(self.gene_ps_results, cluster + ending)):
                pass
            else:
                error_list.append("[!]\t{} does not exist".format(cluster + ending))
        return error_list

    def load_data_and_initialize_global_variables(self):
        def mod_next():
            return next(mg).strip().split(":")
        cluster_count = 0
        input_scope = 0
        error_list = []
        for subdir, dirs, files in os.walk(self.gene_ps_results):
            for file_path_str in files:
                if file_path_str.split(".")[-1] == "GenePS" and len(file_path_str.split(".")) > 1:
                    group_file = os.path.join(subdir, file_path_str)
                    with open(group_file) as mg:
                        group_name = mg.readline().split(":")[1].strip()
                        self.group_names.append(group_name)
                        self.group_to_consensus_file[group_name] = os.path.join(self.gene_ps_results, group_name + ".fa.consensus")
                        self.group_to_original_path[group_name] = group_file
                        self.group_to_group_size[group_name] = int(mg.readline().split(":")[1].strip())
                        input_scope += self.group_to_group_size[group_name]
                        for line in mg:
                            if line.startswith("#name:"):
                                cluster = line.split(":")[1].strip()
                                cluster_count += 1
                                score = mod_next()[1].strip().split(",")
                                self.group_by_cluster_to_score_cutoff[group_name][cluster] = float(score[0])
                                self.group_by_cluster_to_length_range[group_name][cluster] = mod_next()[1].strip().split(",")
                                if os.path.exists(os.path.join(self.gene_ps_results, cluster + ".TN_hmm")):
                                    self.group_by_cluster_to_TN_hmm[group_name][cluster] = os.path.join(self.gene_ps_results, cluster + ".TN_hmm")
                                else:
                                    self.group_by_cluster_to_TN_hmm[group_name][cluster] = None
                                files_not_found = self.validate_path_files(cluster)
                                if not files_not_found:
                                    self.group_by_cluster_to_hmm[group_name][cluster] = os.path.join(self.gene_ps_results, cluster + ".hmm")
                                    self.group_by_cluster_to_fasta_file[group_name][cluster] = os.path.join(self.gene_ps_results, cluster + ".fasta")
                                    self.group_by_cluster_to_fasta_hash[group_name][cluster] = hash_fasta(self.group_by_cluster_to_fasta_file[group_name][cluster])
                                else:
                                    error_list += files_not_found
                            else:
                                pass
        return self.check_loaded_data(cluster_count, input_scope, error_list)


########################################################################################################################
# Class Overseer
########################################################################################################################

class Overseer:
    def __init__(self, g_prefix, genome_location, prediction_location):
        self.g_prefix = g_prefix
        self.genome_path = genome_location
        self.root_directory = prediction_location
        self.group_to_out_dir = {}
        self.group_to_blast_obj = {}
        self.exonerate_file_paths = []

        # statistics
        self.input_scope = 0
        self.merged_regions = 0
        self.filter_count = 0
        self.group_to_fragmented_predictions = defaultdict(int)
        self.group_by_cluster_to_amount_correctly_predicted = defaultdict(lambda: defaultdict(int))

        # results
        self.group_by_contig_to_passed_prediction_list = defaultdict(lambda: defaultdict(list))
        self.group_by_cluster_by_contig_to_filtered_prediction = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        self.group_by_cluster_by_contig_to_valid_prediction = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    ####################################################################################################################
    # Overseer: Initialization
    ####################################################################################################################

    def make_group_directories(self):
        for group in data_base.group_names:
            if self.g_prefix != "single_genome":
                folder_name = self.g_prefix + "_" + group
            else:
                folder_name = group
            group_specific_dir = os.path.join(self.root_directory, folder_name)
            if not os.path.exists(group_specific_dir):
                os.mkdir(group_specific_dir)
            else:
                logging.warning("[!] {} folder already exists; content will be overwritten".format(folder_name))
            self.group_to_out_dir[group] = os.path.join(group_specific_dir, folder_name)
        return self.group_to_out_dir

    ####################################################################################################################
    # Overseer: Find Region
    ####################################################################################################################

    # --> from_header_list_to_fasta
    def blast_all_consensus(self, tmp_directory):
        self.merged_regions = 0
        for group in data_base.group_names:
            consensus_file = data_base.group_to_consensus_file[group]
            if keep:
                blast_obj = run_tblastn(self.genome_path, consensus_file, self.group_to_out_dir[group] + "_intermediate_blast.txt")
            else:
                blast_obj = run_tblastn(self.genome_path, consensus_file, os.path.join(tmp_directory, group))
            if blast_obj is not None:
                blast_obj.infer_regions()
                self.merged_regions += blast_obj.amount_regions
                self.group_to_blast_obj[group] = blast_obj
                if keep:
                    with open(self.group_to_out_dir[group] + "_intermediate_merged_blast_regions.txt", "w") as merged_regions:
                        merged_regions.write(write_merged_region_to_intermediate(blast_obj))
            else:
                logger_blast_region.warning("No Candidate regions found - group: {}".format(group))
        return self.merged_regions

    ##############################################################################################
    # run exonerate on all region -> regional prediction filter -> contigwise_overplapping_control
    ##############################################################################################

    def inform_overseer_about_status(self, status, group, pred_obj):
        if status is True:
            self.group_by_contig_to_passed_prediction_list[group][pred_obj.contig].append(pred_obj)
            return 0
        elif status is False:
            self.group_by_cluster_by_contig_to_filtered_prediction[group][pred_obj.cluster][pred_obj.contig].append(pred_obj)
            logger_validate.info("Filtered - genome: {} group: {} cluster: {} loci: {} {}_{} {} score: {}".format(self.g_prefix, group, pred_obj.cluster, pred_obj.contig, pred_obj.gene_start, pred_obj.gene_end, pred_obj.strand, pred_obj.score))
            return 1
        elif status == "FRAG":
            self.group_to_fragmented_predictions[group] += 1
            self.group_by_cluster_by_contig_to_filtered_prediction[group][pred_obj.cluster][pred_obj.contig].append(pred_obj)
            logger_validate.info("Fragmented - genome: {} group: {} cluster: {} loci: {} {}_{} {} score: {}".format(self.g_prefix, group, pred_obj.cluster, pred_obj.contig, pred_obj.gene_start, pred_obj.gene_end, pred_obj.strand, pred_obj.score))
            return 1
        else:
            raise SyntaxError("[!] only TRUE FALSE FRAG possible in 'inform_overseer_about_status")

    def inform_overseer_about_overlap(self, group, filtered_obj_list, passed_obj_list):
        filtered = 0
        for f_pred in filtered_obj_list:
            filtered += 1
            self.group_by_cluster_by_contig_to_filtered_prediction[group][f_pred.cluster][f_pred.contig].append(f_pred)
            logger_validate.info("overlapping predictions - genome: {} group: {} cluster: {} loci: {} {}_{} {}".format(self.g_prefix, group, f_pred.cluster, f_pred.contig, str(f_pred.gene_start), str(f_pred.gene_end), f_pred.strand))
        for p_pred in passed_obj_list:
            print("\t[+] VALID PREDICTION - genome: {} group: {} cluster: {} loci: {} {}_{} {} score: {}".format(self.g_prefix, group, p_pred.cluster, p_pred.contig, p_pred.gene_start, p_pred.gene_end, p_pred.strand, p_pred.score))
            self.group_by_cluster_to_amount_correctly_predicted[group][p_pred.cluster] += 1
            self.group_by_cluster_by_contig_to_valid_prediction[group][p_pred.cluster][p_pred.contig].append(p_pred)
        return filtered

    def get_exonerate_models(self, out_directory):
        for group in self.group_to_blast_obj:
            for contig in self.group_to_blast_obj[group].inferred_regions:
                for cluster, region_list in self.group_to_blast_obj[group].inferred_regions[contig].items():
                    for region in region_list:
                        if coverage_filter(region) is True:
                            region_fasta = self.group_to_blast_obj[group].region_tuple_to_fasta[region]
                            exo_obj, pred_obj_list = find_best_exonerate_result(region, region_fasta, group, cluster, out_directory)
                            if exo_obj is None:
                                self.filter_count += 1
                                logger_prediction.info("No Exonerate prediction - genome: {} group: {} cluster: {} loci: {} {}_{}".format(self.g_prefix, group, cluster, region.contig, region.s_start, region.s_end))
                            else:
                                self.merged_regions += len(pred_obj_list) - 1
                                self.exonerate_file_paths.append(exo_obj.path)
                                for pred_obj in pred_obj_list:
                                    status = prediction_filter(pred_obj, data_base.group_by_cluster_to_TN_hmm[group][cluster])
                                    self.filter_count += self.inform_overseer_about_status(status, group, pred_obj)
                        else:
                            self.filter_count += 1
                            logger_validate.info("Low Coverage - genome: {} group: {} cluster: {} loci: {} {}_{} {} cov_chunk: {} cov_total:{}".format(self.g_prefix, group, cluster, contig, region.s_start, region.s_end, region.strand, region.chunk_cov, region.query_cov))
                overlap_list, passed_list = isolate_overlapping_predictions(self.group_by_contig_to_passed_prediction_list[group][contig])
                self.filter_count += self.inform_overseer_about_overlap(group, overlap_list, passed_list)
        return self.merged_regions - self.filter_count

    ####################################################################################################################
    # Overseer: write to files
    ####################################################################################################################

    #############################
    # writing fasta and gff files
    #############################

    def write_fasta_gff_files(self, group, group_cluster_contig_prediction, output_type="VALID"):
        protein_dna_gff_array = [[], [], []]
        file_name_list = ["_{}_protein.fa", "_{}_dna.fa", "_{}.gff"]
        file_count = 0
        for cluster in group_cluster_contig_prediction[group]:
            for contig in group_cluster_contig_prediction[group][cluster]:
                for p_obj in group_cluster_contig_prediction[group][cluster][contig]:
                    fasta_header = ">{} geneID:{} Location:{};{}-{} HMM_score:{} Strand:{}".format(
                        genome_prefix, p_obj.geneID, contig, str(p_obj.gene_start), str(p_obj.gene_end), str(round(p_obj.score, 3)), p_obj.strand)
                    protein_dna_gff_array[0].extend([fasta_header, p_obj.protein])
                    protein_dna_gff_array[1].extend([fasta_header, p_obj.DNA])
                    protein_dna_gff_array[2].append("\n".join(p_obj.gff))
                    file_count += 1
        for idx in range(0, len(file_name_list)):
            with open(self.group_to_out_dir[group] + file_name_list[idx].format(output_type), "w") as out_f:
                out_f.write("\n".join(protein_dna_gff_array[idx]))
        return file_count

    #########################################################
    # writing non statistic output --> write_to_seq_gff_files
    #########################################################

    def write_output(self):
        written_valid_files = 0
        written_filtered_files = 0
        for group in self.group_to_out_dir:
            written_valid_files += self.write_fasta_gff_files(group, self.group_by_cluster_by_contig_to_valid_prediction, output_type="VALID")
            written_filtered_files += self.write_fasta_gff_files(group, self.group_by_cluster_by_contig_to_filtered_prediction, output_type="FILTERED")
            if keep:
                with open(self.group_to_out_dir[group] + "_intermediate_exonerate.txt", "wb") as exo_file:
                    for exof in self.exonerate_file_paths:
                        with open(exof, "rb") as infile:
                            exo_file.write(infile.read())
        return written_valid_files, written_filtered_files

    ########################
    # get statistics summary
    ########################

    def get_duplicated_cluster_amount(self, group):
        cluster_count_dict = self.group_by_cluster_to_amount_correctly_predicted
        return len([cluster_count_dict[group][cluster] for cluster in cluster_count_dict[group] if cluster_count_dict[group][cluster] > 1])

    def get_summary_statistics(self):
        group_array = [["group", "#proteins", "#found", "%missing", "#fragmented", "#>1"]]
        cluster_array = [["group", "protein", "contig", "start", "end", "strand", "score"]]
        for group in data_base.group_names:
            if group in self.group_by_cluster_by_contig_to_valid_prediction:
                found = 0
                for cluster in data_base.group_by_cluster_to_score_cutoff[group]:
                    if cluster not in self.group_by_cluster_by_contig_to_valid_prediction[group]:
                        cluster_array.append([group, cluster] + ["-"] * 5)
                    else:
                        found += 1
                        for contig, pred_obj_list in self.group_by_cluster_by_contig_to_valid_prediction[group][cluster].items():
                            for pred_obj in pred_obj_list:
                                cluster_array.append([group, cluster, contig, pred_obj.gene_start, pred_obj.gene_end, pred_obj.strand, round(pred_obj.score, 3)])
                fragmented = self.group_to_fragmented_predictions[group]
                missing = round(((data_base.group_to_group_size[group] - found)/data_base.group_to_group_size[group]) * 100, 2)
                group_array.append([group, data_base.group_to_group_size[group], found, missing, fragmented, self.get_duplicated_cluster_amount(group)])
            else:
                group_array.append([group, data_base.group_to_group_size[group], "0", data_base.group_to_group_size[group]] + ["-"] * 2)
        return convert_summary_to_string(group_array), convert_summary_to_string(cluster_array)

########################################################################################################################
# Function to invoke Overseer genome_wise
########################################################################################################################


def run_GenePS_on_single_genome(current_genome, genome_location, mode="exonerate"):
    prediction_output = get_outdir(out_dir, add_dir=current_genome)
    overseer_obj = Overseer(current_genome, genome_location, prediction_output)
    overseer_obj.make_group_directories()
    if mode == "exonerate":
        db_path = make_blast_db(genome_location, os.path.split(genome_location)[0])
        amount_merged_regions = overseer_obj.blast_all_consensus(tmp_dir)
        amount_valid_predictions = overseer_obj.get_exonerate_models(tmp_dir)
    else:
        raise Exception("[!] unknown mode: {}".format(mode))
    written_valid, written_filtered = overseer_obj.write_output()
    if written_valid == amount_valid_predictions and overseer_obj.filter_count == (overseer_obj.merged_regions - amount_valid_predictions):
        return overseer_obj.get_summary_statistics()
    else:
        print(written_filtered, written_valid, overseer_obj.merged_regions, overseer_obj.filter_count, amount_valid_predictions)
        print("\n[!] unexpected proportions\n")
        logging.error("number of Clusters written to files is not in line with the expectations")
        sys.exit()

###############
# main function
###############

if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)

    print("\n[+] Checking Arguments and Dependencies...")
    check_programs("tblastn", "makeblastdb", "exonerate")
    genome_dict = check_arguments(args)
    print("   -> Okay\n")
    print("[+] Reading GenePS-Files...")
    data_base = DataProviderObject(gene_ps_results)
    cluster_scope = data_base.cluster_scope
    print("   -> Okay\n")

    print("\n\n# {} Genome(s) - {} Group(s) - {} Cluster(s)\n{}\n".format(str(len(genome_dict)), str(len(data_base.group_names)), str(cluster_scope), "#"*42))
    log_path = os.path.join(out_dir, "LOG.txt")
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', filename=log_path, filemode='w')
    if verbose is not None:
        formatter = logging.Formatter('\t%(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)
    logging.info("# Input Directory: {} \n\n".format(gene_ps_results))

    ####################
    # run on all genomes
    ####################
    with tempdir() as tmp_dir:
        genome_count = 1
        for genome_path, genome_prefix in genome_dict.items():
            logging.info("[{}] Analysing genome: {}\n".format(str(genome_count), genome_prefix))
            group_summary, cluster_summary = run_GenePS_on_single_genome(genome_prefix, genome_path, mode="exonerate")
            with open(os.path.join(out_dir, "{}/summary.txt".format(genome_prefix)), "w") as sum_file:
                sum_file.write(group_summary)
                sum_file.write("\n{}\n".format(100 * "-"))
                sum_file.write(cluster_summary)
            genome_count += 1
