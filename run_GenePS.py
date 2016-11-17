#!/usr/bin/env python3

'''
Usage: run_GenePS.py                          -m <DIR> -g <FILE> [-c <INT>] [-s <INT>] [-o <DIR>] [--keep] [--verbose]

    Options:
        -h, --help                            show this screen.

        General
        -m, --GenePS_result_dir <DIR>         folder with consensus/score results from make_GenePS.py
        -g, --genome <FILE|LIST>              Target genome or a ".genomes" file in style of: path/to/genome <TAB> name_prefix
        Optional
        -c, --coverage_filer <INT>            Minimal aligned length of a Blast query to the target genome (used to filter Blast hits)
        -o, --out_dir <DIR>                   Directory for the output files and folders (default: same as input directory)
        -s, --HMM_filter <INT>                Factor to multiply standard deviation of the HMM score distribution with (validate predictions)
        --keep                                Will store intermediate file (Blast output, merged regions, exonerate output)
        --verbose                             Prints log-file information to the screen

'''

import os
import sys
import logging
import tempfile as tmp
from docopt import docopt
from run_command import tempdir, check_programs
from exonerate_parser import run_exonerate, grap_values
from find_regions import run_tblastn, make_blast_db
from collections import defaultdict
from make_GenePS import get_phmm_score, write_to_tempfile, get_outdir

########################################################################################################################
# Global Functions
########################################################################################################################
std_factor = 2
coverage_min = 30
data_base = None
gene_ps_results = None
out_dir = None
keep = None
verbose = None
genome = None
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logger_blast_region = logging.getLogger("BLAST")
logger_prediction = logging.getLogger("Exonerate")
logger_validate = logging.getLogger("Filtering")


def check_arguments(args):
    global coverage_min, std_factor, out_dir, gene_ps_results, keep, verbose, genome
    gene_ps_results = os.path.abspath(args['--GenePS_result_dir'])
    genome_hash = None
    out_dir = "/".join(gene_ps_results.split("/")[:-1])
    if args['--keep']:
        keep = args['--keep']
    if args['--verbose']:
        verbose = args['--verbose']
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
    if args['--HMM_filter']:
        try:
            std_factor = float(args['--HMM_filter'])
        except ValueError:
            error_list.append("[!]\t ERROR: coverage_min needs float or integer; '{}'".format(std_factor))
    if args["--out_dir"]:
        if os.path.isdir(args["--out_dir"]):
            out_dir = os.path.abspath(args["--out_dir"])
        else:
            error_list.append("[!]\t ERROR: output directory: {} does not exist".format(out_dir))
    if args['--genome'].split(".")[-1] == "txt":
        genome_hash = {}
        with open(args['--genome']) as g_file:
            for genome_line in g_file:
                try:
                    genome_line = genome_line.strip("\n").split("\t")
                    len(genome_line) == 2
                    if os.path.exists(genome_line[0]):
                        genome_hash[genome_line[0]] = genome_line[1]
                    else:
                        error_list.append("[!]\t ERROR: {} does not exist or was specified wrongly".format(genome_line[0]))
                except IndexError:
                    error_list.append("[!]\t ERROR: genome file is not a TAB-SEPARATED 2-column file".format(genome_line))
    else:
        error_list.append("[!]\t ERROR: genome file is not a txt file".format(args['--genome']))
    if error_list:
        print("[!] {} - ARGUMENT ERRORS\n".format(len(error_list)))
        print("\n".join(error_list))
        sys.exit()
    return genome_hash


def coverage_filter(area):
    if area.chunk_cov == 100 and area.query_cov < 40:
        return False
    elif area.chunk_cov > coverage_min:
        return True


def convert_summary_to_string(summary_array):
    summary_list = []
    row_format = "\t{}" * (len(summary_array[0]) + 1)     # length of header
    for entry in summary_array:
        summary_list.append(row_format.format("", *entry))
    return "\n".join(summary_list)


def clusterwise_prediction_scoring(fasta_string, length_dict, hmm):
    with tmp.NamedTemporaryFile() as ex_file:
        write_to_tempfile(ex_file.name, fasta_string)
        final_score_hash = get_phmm_score(hmm, ex_file.name, cluster_dict=length_dict)
    if final_score_hash is not None:
        return final_score_hash
    else:
        return None


def make_prediction(q_name, cons_file, dir_path, area, region_seq_string, current_genome, current_group):
    with tmp.NamedTemporaryFile() as reg_file:
        write_to_tempfile(reg_file.name, region_seq_string)
        ex_name = "{}_{}_{}_{}_{}.exon".format(current_genome, current_group, q_name, area.contig, area.s_start, area.s_end)
        ex_obj = run_exonerate(ex_name, dir_path, reg_file.name, cons_file)
        return ex_obj


def write_merged_region_to_intermediate(blast_ob):
    results_list = [
        "# Merging Distance: {}, Flanking Distance {}".format(blast_ob.merging_distance, blast_ob.flanking_distance),
        "# Fields: contig, subject, subject_start, subject_end, strand, chunk_coverage, total_coverage, query_length"]
    for contig in blast_ob.inferred_regions:
        for query in blast_ob.inferred_regions[contig]:
            for region in blast_ob.inferred_regions[contig][query]:
                results_list.append("\t".join([region.contig, query, str(region.s_start), str(region.s_end), str(region.strand), str(region.chunk_cov), str(region.query_cov), str(region.q_len)]))
    return "\n".join(results_list)


########################################################################################################################
# Class: DataProvider (reads in and processes Gene PS files
########################################################################################################################

class DataProviderObject:

    def __init__(self, input_data):
        self.gene_ps_results = input_data
        self.group_names = []
        self.group_to_original_path = {}
        self.group_to_single_copy_ortholog = {}
        self.group_to_group_size = {}
        self.group_by_cluster_to_consensus_length = defaultdict(dict)
        self.group_by_cluster_to_consensus_sequence = defaultdict(dict)
        self.group_by_cluster_to_hmm = defaultdict(dict)
        self.group_by_cluster_to_TN_hmm = defaultdict(dict)
        self.group_by_cluster_to_fasta = defaultdict(dict)
        self.group_by_cluster_to_len_confidence = defaultdict(dict)
        self.group_by_cluster_to_score_cutoff = defaultdict(dict)
        self.cluster_scope = self.load_data_and_initialize_global_variables()

    def check_loaded_data(self, cluster_count, input_scope, error_list):
        if cluster_count != input_scope:
            error_list.append("\t[!] Supposed number of cluster != counted amount of cluster\n")
        if cluster_count == 0:
            error_list.append("\t[!] No GenePS files found - please check your input directory in terms of path and file style\n")
        if error_list:
            print("\n".join(error_list))
            sys.exit()
        else:
            return input_scope

    def validate_path_files(self, cluster):
        error_list = []
        for ending in [".hmmGenePS", ".TN.hmmGenePS", ".fa_GenePS"]:
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
                        self.group_to_original_path[group_name] = group_file
                        self.group_to_group_size[group_name] = int(mg.readline().split(":")[1].strip())
                        input_scope += self.group_to_group_size[group_name]
                        self.group_to_single_copy_ortholog[group_name] = mg.readline().split(":")[1].strip()
                        for line in mg:
                            if line.startswith("#name:"):
                                cluster = line.split(":")[1].strip()
                                cluster_count += 1
                                score = mod_next()[1].strip().split(",")
                                if len(score) > 1:
                                    score = float(score[0]) - (std_factor * float(score[1]))
                                    self.group_by_cluster_to_score_cutoff[group_name][cluster] = score
                                else:
                                    self.group_by_cluster_to_score_cutoff[group_name][cluster] = float(score[0])
                                self.group_by_cluster_to_len_confidence[group_name][cluster] = mod_next()[1].strip().split(",")
                                self.group_by_cluster_to_consensus_sequence[group_name][cluster] = mod_next()[0]
                                self.group_by_cluster_to_consensus_length[group_name][cluster] = (len(self.group_by_cluster_to_consensus_sequence[group_name][cluster]))
                                files_not_found = self.validate_path_files(cluster)
                                if not files_not_found:
                                    self.group_by_cluster_to_hmm[group_name][cluster] = os.path.join(self.gene_ps_results, cluster + ".hmmGenePS")
                                    self.group_by_cluster_to_TN_hmm[group_name][cluster] = os.path.join(self.gene_ps_results, cluster + ".TN.hmmGenePS")
                                    self.group_by_cluster_to_fasta[group_name][cluster] = os.path.join(self.gene_ps_results, cluster + ".fa_GenePS")
                                else:
                                    error_list += files_not_found
                            else:
                                pass
        return self.check_loaded_data(cluster_count, input_scope, error_list)



########################################################################################################################
# Class PredictionObject
########################################################################################################################

class PredictionObject:
    def __init__(self, region, exonerate_obj, score, cutoff, cluster, length_confidence):
        self.cluster = cluster
        self.score = score
        self.cutoff = cutoff
        self.region = region
        self.exon_obj = exonerate_obj
        self.length_confidence = length_confidence
        self.DNA = None
        self.protein = None
        self.gff = None
        self.strand = None
        self.contig = None
        self.blast_location = None
        self.exonerate_location = None
        self.gene_start = None
        self.gene_end = None
        self.gene_length = None
        self.fragmented = None

    def set_up_region_values(self):
        self.DNA = grap_values(self.exon_obj.target_dna)[0]
        self.protein = grap_values(self.exon_obj.target_prot)[0]
        self.gff = grap_values(self.exon_obj.gff)[0]
        self.contig = self.region.contig
        self.strand = self.region.strand
        self.blast_location = (int(self.region.s_start), int(self.region.s_end))
        for target in self.exon_obj.target_prot:
            self.exonerate_location = tuple([int(x) for x in list(self.exon_obj.target_prot[target].keys())[0]])
        self.gene_length = max(self.exonerate_location) - min(self.exonerate_location)
        if self.strand == "-":
            self.gene_start = max(self.blast_location) - min(self.exonerate_location)
            self.gene_end = max(self.blast_location) - max(self.exonerate_location)
        elif self.strand == "+":
            self.gene_start = min(self.blast_location) + min(self.exonerate_location)
            self.gene_end = min(self.blast_location) + max(self.exonerate_location)
        else:
            print("[!]\t ERROR: {} can not be identified as + or - strand".format(self.strand))
            sys.exit()

    def check_for_fragmentation(self):
        if float(self.length_confidence[0]) <= len(self.protein) <= float(self.length_confidence[1]):
            self.fragmented = False
        else:
            self.fragmented = True
        return self.fragmented

    def check_for_overlapp_if_same_contig(self, other):
        if self != other:
            coordinates = [other.gene_start, other.gene_end, self.gene_start, self.gene_end]
            if (other.gene_length + self.gene_length) >= (max(coordinates) - min(coordinates)):
                if self.strand == other.strand:
                    return True
                else:
                    return False
        return False


########################################################################################################################
# Class Overseer
########################################################################################################################

class Overseer:
    def __init__(self, g_prefix, blast_db, prediction_location):
        self.g_prefix = g_prefix
        self.blast_db_path = blast_db
        self.prediction_dir = prediction_location
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
        self.group_by_cluster_blast_filtered_region = defaultdict(lambda: defaultdict(list))
        self.group_by_cluster_no_prediction_region = defaultdict(lambda: defaultdict(list))
        self.group_by_contig_by_position_to_passed_prediction = defaultdict(lambda: defaultdict(list))
        self.group_by_cluster_by_contig_to_filtered_prediction = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        self.group_by_cluster_by_contig_to_valid_prediction = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    ####################################################################################################################
    # Overseer: Initialization
    ####################################################################################################################

    def make_group_directories(self):
        for group in data_base.group_names:
            if self.g_prefix != "single_genome":
                group_name = self.g_prefix + "_" + group
            else:
                group_name = group
            group_specific_dir = os.path.join(self.prediction_dir, group_name)
            if not os.path.exists(group_specific_dir):
                os.mkdir(group_specific_dir)
            else:
                logging.warning("[!] {} folder already existed and content will be overwritten".format(group_name))
            self.group_to_out_dir[group] = os.path.join(group_specific_dir, group_name)
        return self.group_to_out_dir

    ####################################################################################################################
    # Overseer: Find Region
    ####################################################################################################################

    def from_header_list_to_fasta(self, header_list, group_name, file_path):
        with open(file_path + ".consensus", "w") as c_file:
            for header_key in header_list:
                fasta_str = ">{}\n{}\n".format(header_key, data_base.group_by_cluster_to_consensus_sequence[group_name][header_key])
                c_file.write(fasta_str)
        return file_path + ".consensus"

    # --> from_header_list_to_fasta
    def blast_all_consensus(self, out_directory):
        self.merged_regions = 0
        for group in data_base.group_names:
            header_list = list(data_base.group_by_cluster_to_consensus_sequence[group].keys())
            file_path = os.path.join(out_directory, group)
            consensus_file = self.from_header_list_to_fasta(header_list, group, file_path)
            if keep:
                blast_obj = run_tblastn(self.blast_db_path, consensus_file, self.group_to_out_dir[group] + "_intermediate_blast.txt")
            else:
                blast_obj = run_tblastn(self.blast_db_path, consensus_file, os.path.join(out_directory, group))
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

    ####################################################################################################################
    # Overseer: Prediction
    ####################################################################################################################

    ####################
    # make PredictionObj
    ####################

    def clusterwise_benchmarking(self, score_hash, results_hash, cutoff, cluster, group, length_confidence):
        valid_regions = 0
        if score_hash is not None:
            for count, score in score_hash.items():
                prediction_obj = PredictionObject(results_hash[count][1], results_hash[count][0], score, cutoff, cluster, length_confidence)
                prediction_obj.set_up_region_values()
                contig = prediction_obj.contig
                if score >= cutoff:
                    if prediction_obj.check_for_fragmentation() is False:
                        self.group_by_contig_by_position_to_passed_prediction[group][contig].append(prediction_obj)
                        valid_regions += 1
                    else:
                        self.filter_count += 1
                        self.group_to_fragmented_predictions[group] += 1
                        self.group_by_cluster_by_contig_to_filtered_prediction[group][cluster][contig].append(prediction_obj)
                        logger_validate.info("Fragmented Prediction - genome: {} group: {} cluster: {} contig: {} location: {}_{} strand: {}".format(
                            self.g_prefix, group, cluster, contig, prediction_obj.gene_start, prediction_obj.gene_end, prediction_obj.strand))
                else:
                    self.filter_count += 1
                    self.group_by_cluster_by_contig_to_filtered_prediction[group][cluster][contig].append(prediction_obj)
                    logger_validate.info("Filtered Prediction - genome: {} group: {} cluster: {} contig: {} location: {}_{} strand: {} score: {}".format(
                        self.g_prefix, group, cluster, contig, prediction_obj.gene_start, prediction_obj.gene_end, prediction_obj.strand, prediction_obj.score))
        return valid_regions

    #############################################
    # check whether two cluster claim same region
    #############################################

    def contigwise_overlapping_control(self, group, contig):
        count_overlap = 0
        previous_interactions = set()
        for pred_obj in self.group_by_contig_by_position_to_passed_prediction[group][contig]:
            region_owner = pred_obj
            if pred_obj not in previous_interactions:
                highest_score = 0
                for next_pred_obj in self.group_by_contig_by_position_to_passed_prediction[group][contig]:
                    if pred_obj.check_for_overlapp_if_same_contig(next_pred_obj):
                        previous_interactions.update([pred_obj, next_pred_obj])
                        # can still be same cluster but overlapping on different strands!
                        highest_score = max([pred_obj.score, next_pred_obj.score, highest_score])
                        count_overlap += 1
                        if next_pred_obj.score == highest_score:
                            region_owner = next_pred_obj
                            logger_validate.info("overlapping predictions - genome: {} group: {} cluster: {} contig: {} location: {}_{} strand: {}".format(
                                self.g_prefix, group, pred_obj.cluster, contig, str(pred_obj.gene_start), str(pred_obj.gene_end), pred_obj.strand))
                            self.group_by_cluster_by_contig_to_filtered_prediction[group][pred_obj.cluster][contig].append(pred_obj)
                        elif pred_obj.score == highest_score:
                            region_owner = pred_obj
                            self.group_by_cluster_by_contig_to_filtered_prediction[group][next_pred_obj.cluster][contig].append(next_pred_obj)
                            logger_validate.info("overlapping predictions - genome: {} group: {} cluster: {} contig: {} location: {}_{} strand: {}".format(
                                self.g_prefix, group, next_pred_obj.cluster, contig, str(next_pred_obj.gene_start), str(next_pred_obj.gene_end), next_pred_obj.strand))
                        else:
                            pass
                self.group_by_cluster_by_contig_to_valid_prediction[group][region_owner.cluster][region_owner.contig].append(region_owner)
                self.group_by_cluster_to_amount_correctly_predicted[group][region_owner.cluster] += 1
                print("\t[+] VALID PREDICTION - genome: {} group: {} cluster: {} contig: {} location: {}_{} strand: {}".format(
                    self.g_prefix, group, region_owner.cluster, region_owner.contig, region_owner.gene_start, region_owner.gene_end, region_owner.strand))
        return count_overlap

    ############################################################################################
    # run exonerate on all region -> clusterwise_benchmarking -> contigwise_overplapping_control
    ############################################################################################

    def predict_on_all_regions(self, out_directory):
        for group in self.group_to_blast_obj:
            for contig in self.group_to_blast_obj[group].inferred_regions:
                for cluster in self.group_to_blast_obj[group].inferred_regions[contig]:
                    region_count = 1
                    fasta_list = []
                    count_to_exo_to_region, length_hash = {}, {}
                    file_path = os.path.join(out_directory, cluster)
                    consensus_file = self.from_header_list_to_fasta([cluster], group, file_path)
                    len_confidence = data_base.group_by_cluster_to_len_confidence[group][cluster]
                    cut_off = data_base.group_by_cluster_to_score_cutoff[group][cluster]
                    hmm = data_base.group_by_cluster_to_hmm[group][cluster]
                    for region in self.group_to_blast_obj[group].inferred_regions[contig][cluster]:
                        if coverage_filter(region) is True:
                            region_seq = self.group_to_blast_obj[group].region_tuple_to_region_seq[region]
                            exo_obj = make_prediction(cluster, consensus_file, out_directory, region, region_seq, self.g_prefix, group)
                            if exo_obj is None:
                                self.group_by_cluster_no_prediction_region[group][cluster].append((region.contig, region.s_start, region.s_end))
                                self.filter_count += 1
                                logger_prediction.info("No Exonerate prediction - genome: {} group: {} cluster: {} contig: {} region: {}_{}".format(
                                    self.g_prefix, group, cluster, region.contig, region.s_start, region.s_end))
                                continue
                            self.exonerate_file_paths.append(exo_obj.path)
                            predicted_seq = grap_values(exo_obj.target_prot)[0]
                            fasta_header = ">" + str(region_count)
                            count_to_exo_to_region[fasta_header] = (exo_obj, region)
                            length_hash[fasta_header] = len(predicted_seq)
                            fasta_list.extend([fasta_header, predicted_seq])
                        else:
                            self.group_by_cluster_blast_filtered_region[group][cluster].append((region.contig, region.s_start, region.s_end))
                            self.filter_count += 1
                        region_count += 1
                    score_hash = clusterwise_prediction_scoring("\n".join(fasta_list), length_hash, hmm)
                    valid_regions_cluster = self.clusterwise_benchmarking(score_hash, count_to_exo_to_region, cut_off, cluster, group, len_confidence)
                amount_overlapp_contig = self.contigwise_overlapping_control(group, contig)
                self.filter_count += amount_overlapp_contig
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
                    fasta_header = ">{} Cluster:{} Location:{};{}-{} HMM_score:{} Strand:{}".format(
                        genome_prefix, p_obj.cluster, contig, str(p_obj.gene_start), str(p_obj.gene_end), str(round(p_obj.score, 3)), p_obj.strand)
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
                for cluster in data_base.group_by_cluster_to_consensus_length[group]:
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


def locate_all_genes_genome_wise(current_genome, blast_db_path, prediction_out_path):
    overseer_obj = Overseer(current_genome, blast_db_path, prediction_out_path)
    overseer_obj.make_group_directories()
    amount_merged_regions = overseer_obj.blast_all_consensus(tmp_dir)
    amount_valid_predictions = overseer_obj.predict_on_all_regions(tmp_dir)
    written_valid, written_filtered = overseer_obj.write_output()
    if written_valid == amount_valid_predictions and \
                    overseer_obj.filter_count == (amount_merged_regions - amount_valid_predictions):
            return overseer_obj.get_summary_statistics()
    else:
        print("\n[!] unexpected proportions\n")
        logging.error("number of Clusters written to files is not in line with the expectations")
        sys.exit()

###############
# main function
###############

if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)

    #############################
    # read and process input data
    #############################
    print("\n[+] Checking Arguments and Dependencies...")
    check_programs("tblastn", "makeblastdb", "exonerate")
    genome_dict = check_arguments(args)
    print("   -> Okay\n")
    print("[+] Reading GenePS-Files...")
    data_base = DataProviderObject(gene_ps_results)
    cluster_scope = data_base.cluster_scope
    print("   -> Okay\n")

    ###############################
    # make dictionary for all paths
    ###############################
    print("\n[+] Generating BLAST db(s)\n")
    g_prefix_to_db_and_out_dir = {}
    blast_db_dir = get_outdir(out_dir, add_dir="Blast_dbs")
    prediction_dir = get_outdir(out_dir, add_dir="Predictions")
    if genome_dict is None:
        db_path = make_blast_db(genome, blast_db_dir)
        genome_prefix = "single_genome"
        g_prefix_to_db_and_out_dir[genome_prefix] = (db_path, prediction_dir)
    else:
        for genome_path, genome_prefix in genome_dict.items():
            blast_db = get_outdir(blast_db_dir, add_dir=genome_prefix)
            prediction_path = get_outdir(prediction_dir, add_dir=genome_prefix)
            db_path = make_blast_db(genome_path, blast_db)
            g_prefix_to_db_and_out_dir[genome_prefix] = (db_path, prediction_path)

    ######################
    # setting up log files
    ######################
    print("\n\n# {} Genome(s) - {} Group(s) - {} Cluster(s)\n{}\n".format(str(len(g_prefix_to_db_and_out_dir)), str(len(data_base.group_names)), str(cluster_scope), "#"*42))
    log_path = os.path.join(prediction_dir, "LOG.txt")
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
        for single_genome_prefix, db_out_tuple in g_prefix_to_db_and_out_dir.items():
            logging.info("[{}] Analysing genome: {}\n".format(str(genome_count), single_genome_prefix))
            group_summary, cluster_summary = locate_all_genes_genome_wise(single_genome_prefix, db_out_tuple[0], db_out_tuple[1])
            with open(os.path.join(db_out_tuple[1], "summary.txt"), "w") as sum_file:
                sum_file.write(group_summary)
                sum_file.write("\n{}\n".format(100 * "-"))
                sum_file.write(cluster_summary)
            genome_count += 1
