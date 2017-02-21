#!/usr/bin/env python3
import os
from operator import itemgetter
from collections import defaultdict, namedtuple
from shared_code_box import run_cmd


def make_blast_db(genome, temp_dir):
    name = genome.split("/")[-1]
    out_dir = os.path.join(temp_dir, name)
    if os.path.exists(out_dir + ".nhr"):
        print("\n\t[-] BLAST db already exists:\t{}\n".format(os.path.join(temp_dir, name)))
        return out_dir
    command = "makeblastdb -in " + genome + " -dbtype nucl -parse_seqids -out " + out_dir
    run_cmd(command=command, wait=True)
    return out_dir


def parse_blastdb(db_path, contig, start, end):
    region = str(start) + "-" + str(end)
    command = "blastdbcmd -db " + db_path + " -dbtype nucl -entry " + contig + " -range " + region
    results = []
    for line in run_cmd(command=command, wait=False):
        line = line.strip("\n")
        if not line.startswith(">"):
            results.append(line)
    return "".join(results)


def set_min_start(position):
    if position < 0:
        return 1
    else:
        return position


class HspListObject:
    def __init__(self, subject_hsp_list, merging_dist):
        self.merge_dist = merging_dist
        self.hsp_list = subject_hsp_list
        self.hsp_sorted = []
        self.s_start = []
        self.s_end = []
        self.q_len = []
        self.strand = []
        self.q_start = []
        self.q_end = []

    def get_sorted_attributes(self):
        for x in self.hsp_sorted:
            self.s_start.append(x["s_start"])
            self.s_end.append(x["s_end"])
            self.q_len.append(x["q_len"])
            self.strand.append(x["strand"])
            self.q_start.append(x["q_start"])
            self.q_end.append(x["q_end"])

    def sort_hsp_list(self):
        for hsp_x in self.hsp_list:
            if hsp_x["strand"] == "-":
                s_start = hsp_x["s_start"]
                hsp_x["s_start"] = hsp_x["s_end"]
                hsp_x["s_end"] = s_start
        self.hsp_sorted = sorted(self.hsp_list, key=itemgetter("strand", "s_start", "q_start"))
        self.get_sorted_attributes()

    def merge_condition(self, end, start_last, strand, strand_last):
        if abs(end - start_last) < self.merge_dist:
            if strand == strand_last:
                return True
        return False

    def merge_to_region(self):
        if self.hsp_sorted is []:
            self.sort_hsp_list()
        last_idx = 0
        single_merge = None
        all_merged_regions = []
        for idx in range(1, len(self.hsp_sorted)):
            if abs(self.s_end[idx] - self.s_start[last_idx]) < self.merge_dist \
                    and self.strand[idx] == self.strand[last_idx]:
                if single_merge is None:
                    single_merge = [last_idx]
                single_merge.append(idx)
            else:
                if single_merge is not None:
                    all_merged_regions.append(single_merge)
                    single_merge = None
                else:
                    all_merged_regions.append([last_idx])
            last_idx = idx
        if single_merge is not None:
            all_merged_regions.append(single_merge)
        else:
            all_merged_regions.append([len(self.hsp_sorted) - 1])
        return all_merged_regions

    def compute_coverage(self, q_start_pos, q_end_pos, q_length):
        q_starts, q_ends = zip(*sorted(zip(q_start_pos, q_end_pos)))
        aligned_bits_sum = 0
        prev_idx = 0
        for idx in range(1, len(q_ends)):
            if q_starts[idx] < q_ends[prev_idx]:
                if q_ends[idx] < q_ends[prev_idx]:
                    continue
                else:
                    aligned_bits_sum += q_ends[idx] - q_ends[prev_idx]
            else:
                aligned_bits_sum += q_ends[idx] - q_starts[idx]
            prev_idx = idx
        aligned_bits_sum += q_ends[0] - q_starts[0]
        chunck_cov = round((aligned_bits_sum / (max(q_ends) - min(q_starts))) * 100)
        query_cov = round((aligned_bits_sum / q_length) * 100)
        return chunck_cov, query_cov


class BlastObject:
    def __init__(self, blast_out, db_path, merging_dist=10000, flanking_dist=5000):
        self.blast_out = blast_out
        self.db_path = db_path
        self.merging_distance = merging_dist
        self.flanking_distance = flanking_dist
        self.inferred_regions = None
        self.amount_regions = None
        self.region_tuple_to_fasta = {}

    def adjust_oversized_end_pos(self, region_seq, start_position):
        return len(region_seq) + start_position - 1

    def infer_regions(self):
        self.amount_regions = 0
        inferred_regions = {}
        Region = namedtuple('Region', 'contig, s_start, s_end, strand, chunk_cov, query_cov, q_len')
        for subject in self.blast_out:
            if subject not in inferred_regions:
                inferred_regions[subject] = defaultdict(list)
            for query in self.blast_out[subject]:
                hsp_list = self.blast_out[subject][query]
                hits = HspListObject(hsp_list, self.merging_distance)
                hits.sort_hsp_list()
                idx_all_merged_regions = hits.merge_to_region()
                self.amount_regions += len(idx_all_merged_regions)
                for region in idx_all_merged_regions:
                    begin, stop = region[0], region[-1]
                    strand = hits.strand[begin]
                    s_start = set_min_start(hits.s_start[begin] - self.flanking_distance)
                    fictive_s_end = hits.s_end[stop] + self.flanking_distance   # contig may be shorter
                    fasta_str = parse_blastdb(self.db_path, subject, s_start, fictive_s_end)
                    s_end = self.adjust_oversized_end_pos(fasta_str, s_start)
                    q_start_pos = hits.q_start[begin:stop + 1]
                    q_end_pos = hits.q_end[begin:stop + 1]
                    chunk_cov, query_cov = hits.compute_coverage(q_start_pos, q_end_pos, hits.q_len[0])
                    region = Region(contig=subject, s_start=s_start, s_end=s_end, strand=strand,
                                    chunk_cov=chunk_cov, query_cov=query_cov, q_len=hits.q_len[0])
                    inferred_regions[subject][query].append(region)
                    self.region_tuple_to_fasta[region] = ">{}\n{}".format(subject, fasta_str)
        self.inferred_regions = inferred_regions


def read_blast_output(blast_file, db_path):
    blast_dict = {}
    results_flag = False
    with open(blast_file) as blast_f:
        for line in blast_f:
            if not line.startswith("#"):
                results_flag = True
                line = line.strip("\n").split("\t")
                if "-" in line[8]:
                    strand = "-"
                else:
                    strand = "+"
                query = line[0]
                subject = line[1]
                row = {"contig": subject,
                       "evalue": float(line[2]),
                       "q_start": int(line[3]),
                       "q_end": int(line[4]),
                       "s_start": int(line[5]),
                       "s_end": int(line[6]),
                       "q_len": int(line[7]),
                       "strand": strand}
                if subject not in blast_dict:
                    blast_dict[subject] = defaultdict(list)
                blast_dict[subject][query].append(row)
        blast_f.seek(0)
        if results_flag:
            blast_object = BlastObject(blast_dict, db_path)
            return blast_object
        return None


def run_tblastn(db_path, q_file, file_location):
    command = ["tblastn", "-query", q_file, "-db", db_path, "-outfmt",
               """7 qacc sacc evalue qstart qend sstart send qlen sframe""", "-evalue", "1e-1", "-out", file_location]
    run_cmd(command=command, wait=True)
    return read_blast_output(file_location, db_path)

if __name__ == "__main__":
    db = "/home/jgravemeyer/Dropbox/MSc_project/res/c_elegans.PRJNA13758.WS254.genomic.fa"
    db = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/eef_data/F226Dparalog_T15D6.2_region.fa"
    query_test = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/eef_data/25F22D6.12.fa"
    out_dir = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/eef_data/F226Dparalog_T15D6.2.blast"
    blast = run_tblastn(db, query_test, out_dir)
    if blast:
        blast.infer_regions()
    else:
        raise Exception("no results")



    for contig in blast.inferred_regions:
        for query in blast.inferred_regions[contig]:
            print(query)
            for region_x in blast.inferred_regions[contig][query]:
                    print(region_x)
                    print(blast.region_tuple_to_fasta[region_x])

