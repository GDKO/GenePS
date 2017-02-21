#!/usr/bin/env python3
import unittest
import os
from unittest.mock import patch
import build_models
import tempfile as tmp
import shared_code_box
import warnings
warnings.filterwarnings("ignore")

script_path = os.path.dirname(os.path.realpath(__file__))
test_data = os.path.join(script_path, "test_data/compile_script")
file_name = "eef_3.5_no_cea.fa"


class TestHashFasta(unittest.TestCase):

    fa_hash = shared_code_box.hash_fasta(os.path.join(test_data, file_name))

    @patch("build_models.hash_fasta", return_value={">a": 1, ">b": 2})
    @patch("build_models.logger_Filtered", return_value=None)
    def test_check_and_hash_fasta(self, logger, fasta_hash):
        fasta_hash = build_models.check_and_hash_fasta(os.path.join(test_data, file_name), "filename")
        self.assertEqual(fasta_hash, None)

    @patch("build_models.hash_fasta", return_value={})
    @patch("build_models.logger_Filtered", return_value=None)
    def test_check_and_hash_fasta(self, logger, fasta_hash):
        fasta_hash = build_models.check_and_hash_fasta(os.path.join(test_data, file_name), "filename")
        self.assertEqual(fasta_hash, None)

    @patch("build_models.hash_fasta", return_value={">a": 1, ">b": 2, ">c": 3})
    @patch("build_models.logger_Filtered", return_value=None)
    def test_check_and_hash_fasta(self, logger, fasta_hash):
        fasta_hash = build_models.check_and_hash_fasta(os.path.join(test_data, file_name), "filename")
        self.assertEqual(fasta_hash, {">a": 1, ">b": 2, ">c": 3})

    def test_hash_fasta_equal_number_header_seq(self):
        number_seq = len(self.fa_hash.keys())
        number_header = len(self.fa_hash.values())
        self.assertEqual(number_header, number_seq)

    def test_hash_fasta_number_of_entries(self):
        number_entries = len(self.fa_hash.keys())
        self.assertEqual(number_entries, 190)

    def test_hash_fasta_correct_seq_length(self):
        single_seq = len(self.fa_hash[">AMELL.GB42352-PA"])
        self.assertEqual(single_seq, 942)


class TestCleanFastaHash(unittest.TestCase):

    @patch("build_models.logger_Filtered", return_value=None)
    def test_remove_two_keys(self, logger):
        clean_hash = build_models.clean_fasta_hash({">a": 1, ">b": 2, ">c": 3}, [">c"], "test")
        self.assertEqual(clean_hash, {">c": 3})

    @patch("build_models.logger_Filtered", return_value=None)
    def test_key_not_in_hash(self, logger):
        clean_hash = build_models.clean_fasta_hash({">a": 1, ">b": 2, ">c": 3}, [">t", ">p"], "test")
        self.assertEqual(clean_hash, {})

    @patch("build_models.logger_Filtered", return_value=None)
    def test_empty_dictionary(self, logger):
        clean_hash = build_models.clean_fasta_hash({}, [">t", ">p"], "test")
        self.assertEqual(clean_hash, {})


class TestWriteHashToFasta(unittest.TestCase):

    def test_fasta_dictionary(self):
        with tmp.NamedTemporaryFile() as hash_tmp:
            file_path = shared_code_box.write_hash_to_fasta(hash_tmp.name, {">1" : "ATGSAD", ">2": "ADFAT"})
            self.assertSetEqual(set(open(file_path).readlines()), {'>2\n', 'ADFAT\n', '>1\n', 'ATGSAD\n'})

    def test_empty_hash(self):
        with tmp.NamedTemporaryFile() as hash_tmp:
            file_path = shared_code_box.write_hash_to_fasta(hash_tmp.name, {})
            self.assertEqual(file_path, None)

    def test_rigth_optional_formatting(self):
        with tmp.NamedTemporaryFile() as hash_tmp:
            file_path = shared_code_box.write_hash_to_fasta(hash_tmp.name, {">1" : "ATGSAD", ">2": "ADFAT"}, line_style="{}\t{}\t")
            self.assertSetEqual(set(open(file_path).readlines()[0].split("\t")), {'>2', 'ADFAT', '>1', 'ATGSAD', ''})

    def test_formatting_without_dictionary_value(self):
        with tmp.NamedTemporaryFile() as hash_tmp:
            file_path = shared_code_box.write_hash_to_fasta(hash_tmp.name, {">1" : "ATGSAD", ">2": "ADFAT"}, line_style="{}")
            self.assertSetEqual(set(open(file_path).readlines()[0].split(">")), {'', '2', '1'})

    def test_formatting_None_if_more_then_two_brackets(self):
        with tmp.NamedTemporaryFile() as hash_tmp:
            file_path = shared_code_box.write_hash_to_fasta(hash_tmp.name, {">1" : "ATGSAD", ">2": "ADFAT"}, line_style="{}\n{}\n{}\n")
            self.assertEqual(file_path, None)


class TestIntersectionPoint(unittest.TestCase):
    tp = [210.0, 353.0, 403.0, 26.0, 238.0, 409.0, 75.0, 317.0, 279.0, 410.0, 96.0, 311.0, 108.0, 235.0, 375.0, 283.0, 264.0, 193.0, 258.0, 308.0, 355.0, 83.0, 415.0, 22.0, 278.0, 276.0, 442.0, 362.0, 342.0, 450.0, 404.0, 196.0, 372.0, 283.0, 291.0, 194.0, 560.0, 174.0, 251.0, 210.0, 481.0, 272.0, 416.0, 242.0, 313.0, 199.0, 411.0, 322.0, 300.0, 443.0, 237.0, 310.0, 405.0, 223.0, 223.0, 425.0, 227.0, 301.0, 317.0, 166.0, 287.0, 261.0, 253.0, 296.0, 175.0, 311.0, 483.0, 483.0, 141.0, 270.0, 471.0, 390.0, 285.0, 98.0, 279.0, 242.0, 420.0, 119.0, 310.0, 285.0, 412.0, 212.0, 194.0, 279.0, 103.0, 403.0, 228.0, 270.0, 298.0, 283.0, 220.0, 75.0, 411.0, 204.0, 363.0, 283.0, 248.0, 332.0, 326.0, 323.0, 307.0, 95.0, 110.0, 51.0, 336.0, 60.0, 354.0, 248.0, 456.0]
    tn = [257.0, 282.0, 230.0, 147.0, 248.0, 189.0, 248.0, 301.0, 106.0, 271.0, 246.0, 296.0, 138.0, 100.0, 262.0, 166.0]

    def test_normal_intersection(self):
        interpoint = build_models.find_density_intersection(self.tp, self.tn)
        self.assertEqual(interpoint, 300)

    def test_identical_distributions(self):
        interpoint = build_models.find_density_intersection(self.tn, self.tn)
        self.assertEqual(interpoint, None)


class TestLengthBinnedFasta(unittest.TestCase):

    def test_estimate_bin_number_never_below_amount_data(self):
        self.assertEqual(build_models.estimate_bin_number(4, 3, 2, forced_bins=3), 2)

    def test_estimate_bin_force_over_guess(self):
        self.assertEqual(build_models.estimate_bin_number(5, 3, 10, forced_bins=4), 4)

    def test_estimate_bin_emit_guess_if_possible(self):
        self.assertEqual(build_models.estimate_bin_number(5, 3, 10), 5)

    def test_estimate_bin_minimum_before_guess(self):
        self.assertEqual(build_models.estimate_bin_number(5, 6, 10), 6)

    def test_bin_sequence_length_into_three(self):
        bin_dict = build_models.bin_sequence_lengths({"0": 0, "4": 4, "12": 12, "16a": 16, "16b": 16, "18": 18, "24": 24, "26": 26, "28": 28}, 3)
        self.assertEqual(bin_dict, {1: '4', 2: '18', 3: '28'})

    def test_huge_window_size_reduces_bin_number(self):
        self.assertEqual(build_models.bin_sequence_lengths({"0": 0, "4": 4, "12": 12, "16a": 16, "16b": 16, "18": 18, "24": 24, "26": 26, "100": 100}, 3), {1: '26', 3: '100'})

    def test_write_length_binned_fasta_correct_length_hash(self):
        with tmp.NamedTemporaryFile() as tempf:
            length_dict = build_models.write_length_binned_fasta({">first" : "ADFADF", ">second" : "DFAF"}, "test", tempf.name)
            self.assertEqual(length_dict, {'>second': 4, '>first': 6})

    def test_write_length_binned_fasta_converts_from_bin_to_header_into_bin_to_sequence_correctly(self):
        with tmp.NamedTemporaryFile() as tempf:
            length_dict = build_models.write_length_binned_fasta({">first" : "ADFADF", ">second" : "DFAF"}, "test", tempf.name)
            self.assertSetEqual(set(open(tempf.name).readlines()), {'>first_test\n', 'DFAF\n', '>second_test\n', 'ADFADF\n'})


class TestGetpHmmScores(unittest.TestCase):

    hmm_search_score_file = os.path.join(test_data, "domtblout.txt")

    def test_hmmsearch_score_dict_length(self):
        score_dict = shared_code_box.parse_hmmer_domain_table(self.hmm_search_score_file)
        self.assertEqual(len(score_dict), 190)

    def test_hmmsearch_one_domain_score(self):
        score_dict = shared_code_box.parse_hmmer_domain_table(self.hmm_search_score_file)
        self.assertEqual(1583, score_dict[">LOLOA1.EN70_1811"])

    def test_hmmsearch_two_domain_score(self):
        score_dict = shared_code_box.parse_hmmer_domain_table(self.hmm_search_score_file)
        self.assertEqual(1336, score_dict[">SSCAP.L892_g29260.t1"])


class TestScoreObject(unittest.TestCase):

    def setUp(self, size=5):
        self.fasta_dict = {">ACRAS.cds.Contig10403m.5077": "QTFRRIVENINVIIATYGDDDGPMGPIMVDPALGNVGFGSGLHGWAFTLKQFAEMYASKFGVQVDKLMKNLWGDRFFNMKTKKRSTSQEDGAVRGFTQFVLDPIFKVF*",
                  ">ALUMB.ALUE_0000951001-mRNA-1": "XVCVQTETVLRQAIAERIKPVLFMNKMDRALLELQLGQEELYQTFQRIVENTNVIIATYGDDDGPMGQIMVDPAIGNVGFGSGLHGWAFTLKQFAEMYSEKFGVQ",
                  ">ACRAS.cds.Contig3658m.1561": "MEMLYEGPHDDEVAVAIKNCDPNGPLMMYVSKMVPTSDKGRFYAFGRVFSGKVATGMKARIQGPNYVPGKKEDLYEKTIQRTILMMGRYVEPIEDIPSGNIAGLVGVDQYLVKGGTITTFKDAHNLRVMKFSVSPVVRVAVEPKNAGDLPKLVEGLKRLSKSDPMV"}
        self.score_obj = build_models.ScoreObject(self.fasta_dict, os.path.join(test_data, "eef.hmm"))

    def test_iterative_scoring_hand_checked_score(self):
        score_dict = self.score_obj.iterative_score_computation()
        self.assertEqual(score_dict[">ALUMB.ALUE_0000951001-mRNA-1"], 17)
        self.assertNotIn(">ACRAS.cds.Contig3658m.1561", score_dict)

    def test_bulk_scoring_hand_checked_score(self):
        score_dict = self.score_obj.bulk_score_computation()
        self.assertEqual(score_dict[">ALUMB.ALUE_0000951001-mRNA-1"], 190)

    def test_length_range(self):
        length_dict = {">ACRAS.cds.Contig10403m.5077": 10, ">ALUMB.ALUE_0000951001-mRNA-1": 20, ">ACRAS.cds.Contig3658m.1561": 30}
        self.assertListEqual(list(build_models.calculate_length_range(length_dict)), [10.0, 30.0])

    def test_fall_back_to_median_if_weired_distribution(self):
        self.score_obj.score_dict = {x : x for x in range(20, 60, 2)}
        self.assertEqual(20, round(self.score_obj.calculate_score_distribution_parameters(true_negative_scores=list(range(20, 130, 10)))))

    def test_true_intersection_point(self):
        self.score_obj.score_dict = {x:x for x in [210.0, 353.0, 403.0, 26.0, 238.0, 409.0, 75.0, 317.0, 279.0, 410.0, 96.0, 311.0, 108.0, 235.0, 375.0, 283.0, 264.0, 193.0, 258.0, 308.0, 355.0, 83.0, 415.0, 22.0, 278.0, 276.0, 442.0, 362.0, 342.0, 450.0, 404.0, 196.0, 372.0, 283.0, 291.0, 194.0, 560.0, 174.0, 251.0, 210.0, 481.0, 272.0, 416.0, 242.0, 313.0, 199.0, 411.0, 322.0, 300.0, 443.0, 237.0, 310.0, 405.0, 223.0, 223.0, 425.0, 227.0, 301.0, 317.0, 166.0, 287.0, 261.0, 253.0, 296.0, 175.0, 311.0, 483.0, 483.0, 141.0, 270.0, 471.0, 390.0, 285.0, 98.0, 279.0, 242.0, 420.0, 119.0, 310.0, 285.0, 412.0, 212.0, 194.0, 279.0, 103.0, 403.0, 228.0, 270.0, 298.0, 283.0, 220.0, 75.0, 411.0, 204.0, 363.0, 283.0, 248.0, 332.0, 326.0, 323.0, 307.0, 95.0, 110.0, 51.0, 336.0, 60.0, 354.0, 248.0, 456.0]}
        tn = [257.0, 282.0, 230.0, 147.0, 248.0, 189.0, 248.0, 301.0, 106.0, 271.0, 246.0, 296.0, 138.0, 100.0, 262.0, 166.0]
        self.assertEqual(306, self.score_obj.calculate_score_distribution_parameters(true_negative_scores=tn))

    def test_score_distribution_return_half_median(self):
        self.score_obj.score_dict = {x : x for x in range(20, 60, 2)}
        self.assertEqual(20, round(self.score_obj.calculate_score_distribution_parameters()))

    def test_score_distribution_to_less_tp_values(self):
        self.score_obj.score_dict = {x : x for x in range(20, 60, 10)} # just 4
        self.assertEqual(18, round(self.score_obj.calculate_score_distribution_parameters(true_negative_scores=list(range(20, 130, 10)))))

    def test_score_distribution_too_less_tn_values(self):
        self.score_obj.score_dict = {x : x for x in range(20, 60, 2)}
        self.assertEqual(20, round(self.score_obj.calculate_score_distribution_parameters(true_negative_scores=[2, 3, 9])))


class TestMsaObject(unittest.TestCase):

    msa_path = open(os.path.join(test_data, "eef.aln"))
    msa_list = msa_path.readlines()
    msa_path.close()

    @patch("build_models.run_cmd", return_value=msa_list)
    def test_msa_operations(self, run_cmd_mock):
        msa_list = build_models.msa_operations("blaaa")
        self.assertEqual(len(msa_list)/2, 190)


class Overseer(unittest.TestCase):

    @patch("build_models.output_dir", return_value="")
    @patch("builtins.print", autospec=True, return_value=None)
    def test_initalize_with_single_file(self, printfunction, outdir):
        overseer_obj = build_models.Overseer(os.path.join(test_data, file_name))
        overseer_obj.initialize_input_data()
        self.assertEqual(overseer_obj.group_by_file_to_filepath["compile_script"]["eef_3.5_no_cea"], "/home/jgravemeyer/Dropbox/MSc_project/src/GenePS/test_data/compile_script/eef_3.5_no_cea.fa")

    @patch("build_models.output_dir", return_value="")
    @patch("builtins.print", autospec=True, return_value=None)
    @patch("build_models.logger_Filtered", return_value=None)
    def test_initalize_with_directory(self, logger, printfunction, outdir):
        overseer_obj = build_models.Overseer(test_data)
        overseer_obj.initialize_input_data()
        self.assertEqual(overseer_obj.group_by_file_to_filepath["compile_script"]["eef_3.5_no_cea"], "/home/jgravemeyer/Dropbox/MSc_project/src/GenePS/test_data/compile_script/eef_3.5_no_cea.fa")
        self.assertNotIn("eef.hmm", overseer_obj.group_by_file_to_filepath["compile_script"])

    def test_filter_file_list_remove_two(self):
        overseer_obj = build_models.Overseer("test")
        overseer_obj.group_to_file_list["testgroup"] = ["1", "2", "3", "4"]
        overseer_obj.valid_input_scope = 4
        valid_files = overseer_obj.remove_filtered_files({"testgroup" : ["1", "2"]})
        self.assertListEqual(overseer_obj.group_to_file_list["testgroup"], ['3', '4'])
        self.assertEqual(4-2, valid_files)

    def test_filter_file_list_nothing_to_remove(self):
        overseer_obj = build_models.Overseer("test")
        overseer_obj.group_to_file_list["testgroup"] = ["1", "2", "3", "4"]
        overseer_obj.valid_input_scope = 4
        valid_files = overseer_obj.remove_filtered_files({})
        self.assertListEqual(overseer_obj.group_to_file_list["testgroup"], ["1", "2", "3", "4"])
        self.assertEqual(4, valid_files)
'''
    @patch("make_Datasets.output_dir", return_value="")
    @patch("make_Datasets.generate_hmm", return_value="Test_HMM")
    @patch("make_Datasets.write_length_binned_fasta", return_value={"a": 1})
    @patch("make_Datasets.print_progress", return_value="")
    def test_output_hmm_fasta_correct_file_writing_after_filtering(self, progress_print, binned_fasta, hmmbuild, outputdir):
        overseer_obj = make_Datasets.Overseer(os.path.join(test_data, file_name))
        overseer_obj.initialize_input_data()
        print(overseer_obj.group_by_file_to_filepath)
        print(overseer_obj.group_to_file_list)
        print(len(overseer_obj.group_by_file_to_cluster_hash))
        with tempdir() as msa_dir:
            overseer_obj.output_HMM_and_fasta(msa_dir)
            print(overseer_obj.group_by_file_to_msa_obj['compile_script']['eef_3.5_no_cea'].size_history)
        # check if same path holds later on a fasta file and now aln
            # aln afer realining same size as history as well as the length of the purged fasta hash
'''
# 190, 129
'''
def remove_filtered_files(self, removal_group_to_file_list):
    for group, file_list in removal_group_to_file_list.items():
        for file_name in file_list:
            self.group_to_file_list[group].remove(file_name)
            self.valid_input_scope -= 1
    return self.valid_input_scope
'''

if __name__ == '__main__':
    unittest.main()



