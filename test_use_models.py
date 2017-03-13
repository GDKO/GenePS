#!/usr/bin/env python3
import unittest
import os
from unittest.mock import patch
from shared_code_box import tempdir
import use_models
from collections import defaultdict, namedtuple
from Blast_wrapper import read_blast_output, HspListObject
from Exonerate_GenBlast_Wrapper import remove_non_letter_signs, clear_hashed_bases, aacode_3to1, ExonerateObject
import Exonerate_GenBlast_Wrapper

script_path = os.path.dirname(os.path.realpath(__file__))
test_data = os.path.join(script_path, "test_data")


class TestBlastObject(unittest.TestCase):

    blast_file = os.path.join(test_data, "run_geneps/remanei_OG1685.blast")
    blast = read_blast_output(blast_file, os.path.join(test_data, "databases/c_brenneri.PRJNA20035.WS249.genomic.fa"))
    blast.infer_regions()

    def test_inferred_regions_single_hsp_in_hsp_list(self):
        single_hsp_region = self.blast.inferred_regions["Cbre_Contig419"]["OrthologousGroups_I3.5.OGoverlapp.txt"][0]
        self.assertEqual(single_hsp_region.s_start, 60160)
        self.assertEqual(single_hsp_region.s_end, 70339)
        self.assertEqual(single_hsp_region.strand, "+")
        self.assertEqual(single_hsp_region.chunk_cov, 100)
        self.assertEqual(single_hsp_region.query_cov, 11)
        self.assertEqual(single_hsp_region.q_len, 466)

    def test_inferred_regions_merging_reverse_compliment(self):
        single_hsp_region = self.blast.inferred_regions["Cbre_Contig6"]["OrthologousGroups_I3.5.OGoverlapp.txt"][0]
        self.assertEqual(single_hsp_region.s_start, 1479061)
        self.assertEqual(single_hsp_region.s_end, 1492025)
        self.assertEqual(single_hsp_region.strand, "-")
        self.assertEqual(single_hsp_region.chunk_cov, 82)
        self.assertEqual(single_hsp_region.query_cov, 82)
        self.assertEqual(single_hsp_region.q_len, 466)

    def test_inferred_regions_merging_plus_strand(self):
        single_hsp_region = self.blast.inferred_regions["Cbre_Contig231"]["OrthologousGroups_I3.5.OGoverlapp.txt"][0]
        self.assertEqual(single_hsp_region.s_start, 10349)
        self.assertEqual(single_hsp_region.s_end, 24214)
        self.assertEqual(single_hsp_region.strand, "+")
        self.assertEqual(single_hsp_region.chunk_cov, 82)
        self.assertEqual(single_hsp_region.query_cov, 82)
        self.assertEqual(single_hsp_region.q_len, 466)

    def test_unmerged_blast_output(self):
        single_hsp_region = self.blast.blast_out["Cbre_Contig419"]["OrthologousGroups_I3.5.OGoverlapp.txt"][0]
        self.assertEqual(single_hsp_region['s_start'], 65160)
        self.assertEqual(single_hsp_region['s_end'], 65339)

    def test_unmerged_blast_reverse_strand(self):
        blast = read_blast_output(self.blast_file, "fake_db")
        single_hsp_region = blast.blast_out["Cbre_Contig6"]["OrthologousGroups_I3.5.OGoverlapp.txt"][0]
        self.assertLess(single_hsp_region['s_end'], single_hsp_region['s_start'])

    def test_compare_merged_unmerged_number_regions(self):
        blast = read_blast_output(self.blast_file, "fake_db")
        number_unmerged_regions = len(blast.blast_out["Cbre_Contig231"]["OrthologousGroups_I3.5.OGoverlapp.txt"])
        self.assertEqual(number_unmerged_regions, 4)
        number_merged_regions = len(self.blast.inferred_regions["Cbre_Contig231"]["OrthologousGroups_I3.5.OGoverlapp.txt"])
        self.assertEqual(number_merged_regions, 1)

    def test_compute_coverage(self):
        test_obj = HspListObject("fake", "fake")
        chunk, query = test_obj.compute_coverage([35, 40, 50, 80], [50, 60, 55, 100], 100)
        self.assertEqual(chunk, 69)
        self.assertEqual(query, 45)
        chunk_many, query_many = test_obj.compute_coverage([88, 64, 238, 93, 241, 64, 101], [184, 101, 400, 216, 360, 101, 215], 466)
        self.assertEqual(chunk_many, 93)
        self.assertEqual(query_many, 67)
        chunk_single, query_single = test_obj.compute_coverage([35], [50], 100)
        self.assertEqual(chunk_single, 100)
        self.assertEqual(query_single, 15)

'''
class TestExonerateObject(unittest.TestCase):

    # with hashes in sequence
    exonerate_obj = ExonerateObject(os.path.join(test_data, "run_geneps/elegans_eef_true.exonerate"))
    target_hash = "I;9859005-9873867"
    target_range_hash = ('4791', '9857')
    target_prot_hash = exonerate_obj.target_prot[target_hash][target_range_hash]
    target_dna_hash = exonerate_obj.target_dna[target_hash][target_range_hash]

    # no hashes in sequence
    target2 = 'I;9159307-9171689'
    target_range2 = ("819", "10142")
    target_prot2 = exonerate_obj.target_prot[target2][target_range2]
    target_dna2 = exonerate_obj.target_dna[target2][target_range2]

    def test_target_target_range(self):
        self.assertIn(self.target_hash, self.exonerate_obj.query_prot)
        self.assertIn(self.target_range_hash, list(self.exonerate_obj.query_prot[self.target_hash].keys()))

    def test_sequences_dividability(self):
        self.assertEqual(len(self.target_dna_hash) / 3, len(self.target_prot_hash))
        self.assertEqual(len(self.target_dna2) / 3, len(self.target_prot2))

    def test_protein_length(self):
        self.assertEqual(827, len(self.target_prot_hash))
        self.assertEqual(961, len(self.target_prot2))

    def test_dna_length(self):
        self.assertEqual(2481, len(self.target_dna_hash))
        self.assertEqual(2883, len(self.target_dna2))

    def test_gff_before_rewriting(self):
        first_element = 'I;9859005-9873867\texonerate:protein2genome:bestfit\tgene\t4792\t9857\t485\t+\t.\tgene_id 0 ; sequence eef_3.5 ; gene_orientation +'
        self.assertEqual(self.exonerate_obj.gff[self.target_hash][self.target_range_hash][0], first_element)

    def test_gff_cds_phase(self):
        self.assertEqual(self.exonerate_obj.gff[self.target_hash][self.target_range_hash][11].split("\t")[-2], '1')

    def test_gff_mRNA(self):
        self.assertEqual(self.exonerate_obj.gff[self.target_hash][self.target_range_hash][1].split("\t")[2], 'mRNA')

    def test_gff_no_splice_no_similarity(self):
        self.assertNotIn("similarity", self.exonerate_obj.gff[self.target_hash][self.target_range_hash][-1].split("\t")[2])
        self.assertIn("exon", self.exonerate_obj.gff[self.target_hash][self.target_range_hash][-1].split("\t")[2])
        self.assertNotIn("splice", self.exonerate_obj.gff[self.target_hash][self.target_range_hash][3].split("\t")[2])

    def test_aa3_to_1_coding_dict(self):
        test_string = 'Val***XaaUnkSer'
        out_string = 'VXXXS'
        failed_string = "ValUnk*"
        self.assertEqual(out_string, aacode_3to1(test_string))
        self.assertRaises(Exception, aacode_3to1, failed_string)

    def test_clear_hashed_bases(self):
        testprot = "        sHisValGluLeuLeu##LysSer------------------SerGlyIleSerLeuLeuCy"
        testdna = " 5030 : ACACGTGGAATTACTATGAAATCC------------------AGTGGAATCAGTTTGCTCTG : 5073"
        out_dna = " 5030 : ACACGTGGAATTACTAAAATCC------------------AGTGGAATCAGTTTGCTCTG : 5073"
        self.assertEqual(out_dna, clear_hashed_bases(testprot, testdna))

    def test_remove_non_letter_signs(self):
        test_string = "AdfA%y*8l#+-w"
        out_string = "AdfAy*lw"
        self.assertEqual(remove_non_letter_signs(test_string), out_string)

'''


class TestOverlapp(unittest.TestCase):

    def make_pred_object_list(self):
        Region = namedtuple("Region", "contig, s_start, s_end, strand")
        region = Region(contig='II', s_start=3290626, s_end=3305673, strand='+')
        obj_values = [[3300015, 3304612, '+', 44, 104], [3295630, 3297975, '+', 628, 886], [3295630, 3296549, '+', 326, 569],
                      [3298963, 3299551, '+', 194, 362], [3298841, 3299485, '+', 217, 264], [3294179, 3300608, '+', 178, 279],
                      [3297119, 3300674, '+', 579, 790]]
        pred_obj_list = []
        for value_list in obj_values:
            prec_obj = Exonerate_GenBlast_Wrapper.PredictionObject(region, value_list[3], "-", "-", "-")
            prec_obj.gene_start = value_list[0]
            prec_obj.gene_end = value_list[1]
            prec_obj.gene_length = prec_obj.gene_end - prec_obj.gene_start
            prec_obj.aln_score = value_list[-1]
            pred_obj_list.append(prec_obj)
        return sorted(pred_obj_list, key=lambda x: x.score, reverse=True)

    @patch('Exonerate_GenBlast_Wrapper.PredictionObject.check_overlap', return_value=True)
    def test_all_overlapp(self, overlapp):
        filtered, passed = use_models.isolate_overlapping_predictions(self.make_pred_object_list())
        self.assertEqual(passed[0].score, 628)
        self.assertEqual(len(filtered), 6)

    @patch('Exonerate_GenBlast_Wrapper.PredictionObject.check_overlap', return_value=False)
    def test_no_overlapp(self, overlapp):
        filtered, passed = use_models.isolate_overlapping_predictions(self.make_pred_object_list())
        self.assertEqual(len(passed), 7)
        self.assertEqual(len(filtered), 0)

    #@patch('Exonerate_GenBlast_Wrapper.PredictionObject.check_overlap', autospec=True, side_effect=overlapp_side_effect)
    def test_crossoverlapping_correct_pick(self):
        filtered, passed = use_models.isolate_overlapping_predictions(self.make_pred_object_list())
        winner_scores = [628, 217, 44]
        for pred_obj in passed:
            self.assertIn(pred_obj.score, winner_scores)

    def test_crossoverlapping_list_length(self):
        filtered, passed = Exonerate_GenBlast_Wrapper.isolate_overlapping_predictions(self.make_pred_object_list())
        self.assertEqual(len(filtered), 4)
        self.assertEqual(len(passed), 3)
        self.assertEqual(len(set(passed)), 3)

    def test_no_overlapp_between_filtered_and_passed(self):
        filtered, passed = Exonerate_GenBlast_Wrapper.isolate_overlapping_predictions(self.make_pred_object_list())
        self.assertFalse(set(passed).intersection(set(filtered)))

    def test_no_redundance_of_passed_objects(self):
        filtered, passed = Exonerate_GenBlast_Wrapper.isolate_overlapping_predictions(self.make_pred_object_list())
        for idx in range(0, len(passed)):
            for idx3 in range(idx + 1, len(passed)):
                self.assertNotEqual(passed[idx], passed[idx+1])


class TestDataProviderObject(unittest.TestCase):

    data_base = use_models.DataProviderObject(test_data + "/run_geneps")

    def test_all_attributes(self):
        self.assertEqual(self.data_base.cluster_scope, 5)
        self.assertEqual(self.data_base.group_names[0], "next_best_blast_eef")
        self.assertEqual(self.data_base.group_to_group_size["next_best_blast_eef"], 5)
        self.assertTrue(self.data_base.group_by_cluster_to_hmm["next_best_blast_eef"]["eef_3.5"].split(".")[-1] == "hmm")
        self.assertEqual(self.data_base.group_by_cluster_to_score_cutoff["next_best_blast_eef"]["eef_3.5"], 1.1826699418473992)
        self.assertEqual(self.data_base.group_by_cluster_to_length_range["next_best_blast_eef"]["eef_3.5"], ['303.97390927100224', '913.1505146921315'])

'''

class TestOverseerPredictAllRegions(unittest.TestCase):
    elegans_db = os.path.join(test_data, "databases/c_elegans.PRJNA13758.WS254.genomic.fa")
    overseer = run_GenePS.Overseer("c_elegans", elegans_db, "prediction_location")
    blast = read_blast_output(os.path.join(test_data, "run_geneps/elegans_blast_eef_nextBest.txt"), elegans_db)
    blast.infer_regions()
    overseer.group_to_blast_obj["next_best_blast_eef"] = blast
    overseer.merged_regions = 21
    run_GenePS.data_base = run_GenePS.DataProviderObject(test_data + "/run_geneps")
    with tempdir() as tmp:
        valid_predictions = overseer.get_exonerate_models(tmp)

    def test_predict_all_regions_number_valid_predictions(self):
        self.assertEqual(5, self.valid_predictions)
        self.assertEqual(17, self.overseer.filter_count)

    def test_predict_all_regions_correct_cluster(self):
        cluster_list = ["OrthologousGroups_I3.5.OG0000365.txt", "OrthologousGroups_I3.5.OG0002137.txt", "OrthologousGroups_I3.5.OG0000365.txt", "eef_3.5"]
        cluster_test_list = []
        for cluster in self.overseer.group_by_cluster_by_contig_to_valid_prediction["next_best_blast_eef"]:
            cluster_test_list.append(cluster)
        self.assertSetEqual(set(cluster_test_list), set(cluster_list))

    def test_predict_all_regions_minus_strand_prediction(self):
        check_list = ['III', 1193723, 1185388, '-']
        for pred in self.overseer.group_by_cluster_by_contig_to_valid_prediction["next_best_blast_eef"]["OrthologousGroups_I3.5.OG0000365.txt"]["III"]:
            self.assertListEqual(check_list, [pred.contig, pred.gene_start, pred.gene_end, pred.strand])

    def test_predict_all_regions_plus_strand_prediction(self):
        check_list = ['II', 13121212, 13121490, '+']
        for pred in self.overseer.group_by_cluster_by_contig_to_valid_prediction["next_best_blast_eef"]["OrthologousGroups_I3.5.OG0000365.txt"]["II"]:
            self.assertListEqual(check_list, [pred.contig, pred.gene_start, pred.gene_end, pred.strand])


class TestGenBlastObject(unittest.TestCase):

    file_path_dict = {'gff': os.path.join(test_data, "run_geneps/python_genblast_test_1.1c_2.3_s1_0_16_1.gff"),
                      'DNA': os.path.join(test_data, 'run_geneps/python_genblast_test_1.1c_2.3_s1_0_16_1.DNA'),
                      'pro': os.path.join(test_data, 'run_geneps/python_genblast_test_1.1c_2.3_s1_0_16_1.pro')}
    hmm_scores = {'>9_B0205.10_Inf5.0_OG0011445-R1-1-A1': 276, '>17_B0205.10_Inf5.0_OG0011445-R1-1-A1': 283,
                  '>12_B0205.10_Inf5.0_OG0011445-R1-1-A1': 301, '>1_B0205.10_Inf5.0_OG0011445-R1-1-A1': 295,
                  '>16_B0205.10_Inf5.0_OG0011445-R1-1-A1': 259, '>15_B0205.10_Inf5.0_OG0011445-R1-1-A1': 281,
                  '>13_B0205.10_Inf5.0_OG0011445-R1-1-A1': 403, '>20_B0205.10_Inf5.0_OG0011445-R1-1-A1': 292,
                  '>19_B0205.10_Inf5.0_OG0011445-R1-1-A1': 306, '>14_B0205.10_Inf5.0_OG0011445-R1-1-A1': 332,
                  '>8_B0205.10_Inf5.0_OG0011445-R1-1-A1': 320, '>1_B0205.10_Inf5.0_OG0011445-R1-1-A2': 293}

    @patch('Exonerate_GenBlast_Wrapper.run_cmd', return_value=None)
    def test_run_genblastg_returns_3_files(self, run_cmd):
        genblast = Exonerate_GenBlast_Wrapper.run_genblastg("", "", os.path.join(test_data, "run_geneps"), "python_genblast_test")
        self.assertEqual(self.file_path_dict, genblast.out_files)
        pass

    def test_hash_gff_correct_length(self):
        genblast_obj = Exonerate_GenBlast_Wrapper.GenblastObject(self.file_path_dict)
        genblast_obj.hash_gff_and_infer_regions()
        self.assertEqual(13, len(genblast_obj.gff))
        self.assertEqual(13, len(genblast_obj.regions))

    def test_hash_gff_test_phase(self):
        genblast_obj = Exonerate_GenBlast_Wrapper.GenblastObject(self.file_path_dict)
        genblast_obj.hash_gff_and_infer_regions()
        transcript = genblast_obj.gff[">13_B0205.10_Inf5.0_OG0011445-R1-1-A1"][0].split("\t")[7]
        first_exon = genblast_obj.gff[">13_B0205.10_Inf5.0_OG0011445-R1-1-A1"][1].split("\t")[7]
        second_exon = genblast_obj.gff[">13_B0205.10_Inf5.0_OG0011445-R1-1-A1"][2].split("\t")[7]
        third_exon = genblast_obj.gff[">13_B0205.10_Inf5.0_OG0011445-R1-1-A1"][3].split("\t")[7]
        last_exon = genblast_obj.gff[">13_B0205.10_Inf5.0_OG0011445-R1-1-A1"][4].split("\t")[7]
        phase_list = [transcript, first_exon, second_exon, third_exon, last_exon]
        self.assertListEqual(phase_list, [".", "0", "1", "0", "1"])

    def test_region_tupel(self):
        genblast_obj = Exonerate_GenBlast_Wrapper.GenblastObject(self.file_path_dict)
        genblast_obj.hash_gff_and_infer_regions()
        region = genblast_obj.regions[">13_B0205.10_Inf5.0_OG0011445-R1-1-A1"]
        attribute_list = [region.contig, region.s_start, region.s_end, region.strand, region.header]
        self.assertListEqual(attribute_list, ["I", "10708833", "10710342", "+", ">13_B0205.10_Inf5.0_OG0011445-R1-1-A1"])

    def test_infer_prediction_objects(self):
        genblast_obj = Exonerate_GenBlast_Wrapper.GenblastObject(self.file_path_dict)
        genblast_obj.hash_gff_and_infer_regions()
        obj_list = genblast_obj.infer_prediction_objects(self.hmm_scores, "cluster", "cut_off", "length_range")
        for obj in obj_list:
            header = obj.region.header
            self.assertEqual(obj.score, self.hmm_scores[header])
            self.assertEqual(obj.gene_length, int(genblast_obj.regions[header].s_end) - int(genblast_obj.regions[header].s_start))
            self.assertEqual(obj.strand, genblast_obj.regions[header].strand)
            self.assertEqual(obj.contig, genblast_obj.regions[header].contig)
            self.assertEqual(obj.cluster, "cluster")
            self.assertEqual(obj.cutoff, "cut_off")

    def test_remove_overlapping_predictions(self):
        genblast_obj = Exonerate_GenBlast_Wrapper.GenblastObject(self.file_path_dict)
        genblast_obj.hash_gff_and_infer_regions()
        obj_list = genblast_obj.infer_prediction_objects(self.hmm_scores, "cluster", "cut_off", "length_range")
        filtered, passed = genblast_obj.separate_self_overlapping_predictions(self.hmm_scores, obj_list)
        self.assertEqual(len(filtered) + len(passed), len(obj_list))
        for p_obj in passed:
            self.assertIn(p_obj.region.header, [">13_B0205.10_Inf5.0_OG0011445-R1-1-A1", ">14_B0205.10_Inf5.0_OG0011445-R1-1-A1"])

    def test_fill_prediction_objects(self):
        genblast_obj = Exonerate_GenBlast_Wrapper.GenblastObject(self.file_path_dict)
        genblast_obj.hash_gff_and_infer_regions()
        obj_list = genblast_obj.infer_prediction_objects(self.hmm_scores, "cluster", "cut_off", "length_range")
        filtered, passed = genblast_obj.separate_self_overlapping_predictions(self.hmm_scores, obj_list)
        genblast_obj.fill_prediction_objects()
        for test_obj in genblast_obj.pred_obj:
            if test_obj.region.header == ">13_B0205.10_Inf5.0_OG0011445-R1-1-A1":
                self.assertEqual(len(test_obj.protein), 446)
                self.assertEqual(len(test_obj.DNA), 1338)
                self.assertEqual(len(test_obj.gff), 5)
'''

if __name__ == '__main__':
    unittest.main()
