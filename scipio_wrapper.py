#!/usr/bin/env python3
from Exonerate_GenBlast_Wrapper import PredictionObject, isolate_overlapping_predictions
import os
from collections import defaultdict, namedtuple


def get_gff_files(directory):
    file_list = []
    for x,y,z in os.walk(directory):
        for file in z:
            if file.split(".")[-1] == "gff":
                file_list.append(os.path.join(x, file))
    return file_list


def make_gff_dict(gff_file):
    with open(gff_file) as gff:
        count = 0
        filename = ".".join(gff_file.split(".")[:-3]).split("/")[-1]
        gff_dict = defaultdict(list)
        for line in gff:
            if line.startswith("##query"):
                count += 1
            elif line.startswith("#"):
                pass
            else:
                new_line = line.strip("\n").split("\t")
                new_line[2] = "cds"
                new_line[-1] = new_line[-1].split(" ")[0].split("%")[0]
                gff_dict["{}_{}".format(str(count), filename)].append(new_line)
    return gff_dict


def gff_to_pred_obj(gff_dictionary):
    Region = namedtuple('Region', 'strand, contig')
    pred_obj_list = []
    for entry, gff_list in gff_dictionary.items():
        strand = gff_list[0][6]
        contig = gff_list[0][0]
        score = float(gff_list[0][5])
        region = Region(strand=strand, contig=contig)
        pred_obj = PredictionObject(region, score, entry, "cutoff", "lenght_range")
        if strand == "+":
            pred_obj.gene_start = int(gff_list[0][3])
            pred_obj.gene_end = int(gff_list[-1][4])
        else:
            pred_obj.gene_start = int(gff_list[0][4])
            pred_obj.gene_end = int(gff_list[-1][3])
        pred_obj.gene_length = abs(pred_obj.gene_end - pred_obj.gene_start)
        if pred_obj.score >= 0.1:
            pred_obj.gff = ["\t".join(x) for x in gff_list]
            pred_obj_list.append(pred_obj)
    filtered, passed = isolate_overlapping_predictions(pred_obj_list)
    return passed


scipio_dir = "/home/jgravemeyer/programs/scipio-1.4/scipio_accuracy"
output = "/home/jgravemeyer/programs/scipio-1.4/scipio_accuracy/scipio_all_cluster.gff"
scipio_gffs = get_gff_files(scipio_dir)
result_list_gff = []
for gff in scipio_gffs:
    #if "B0205.10_Inf5.0_OG0011445" in gff:
        gff_dict = make_gff_dict(gff)
        result_list_gff.extend(gff_to_pred_obj(gff_dict))
print(len(result_list_gff))

with open(output, "w") as outf:
    outf.write("##ggf-version\t3\n")
    for pred in result_list_gff:
        outf.write("\n".join(pred.gff) + "\n")


### load all gff files
### per entry per gff file predciton object
### per gff file overlap function ( manual assembly yes and no)
### write all gff entries in one new gff file
