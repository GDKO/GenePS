#!/usr/bin/env python3
# plotting true negatives vs true positives

"""
Usage: plot_PHMM_scores.py                      -i <DIR> -o <DIR> [--membership <FILE>]

    Options:
        -h, --help                            show this screen.

        General
        -i, --input_dir <DIR>                 directory with intermediate files from build_models.py
        -o, --out_dir <DIR>                   output directory
        --membership <FILE>                   seeks and plots for corresponding orthofinder groups of next best BLAST proteins
"""

import os
import sys
import numpy as np
from docopt import docopt
from build_models import find_density_intersection
from shared_code_box import run_cmd, get_outdir

import_errors = []
try:
    import matplotlib.pyplot as plt
except ImportError:
    import_errors.append("[ERROR] : Module \'matplotlib.pyplot\' was not found. Please install \'matplotlib\' using your package manager")
try:
    import seaborn as sns
    sns.set(color_codes=True)
except ImportError:
    import_errors.append("[ERROR] : Module \'seaborn\' was not found. Please install \'seaborn\' using \'pip install seaborn\'")
if import_errors:
    sys.exit("\n".join(import_errors))


def grap_output (protein_id, orthologous_file):
    command = ["grep", """{}""".format(protein_id), orthologous_file]
    for line in run_cmd(command=command, wait=False):
        group = line.strip().split(":")[0]
        return group


def hash_score_files(score_file):
    score_hash = {}
    with open(score_file) as sf:
        for line in sf:
            if not line.startswith("#"):
                line = line.strip("\n").split("\t")
                score_hash[line[0]] = round(float(line[1]), 3)
    return score_hash, list(score_hash.values())

'''
def extract_scores(score_dict, normalize=False):
    score_list = []
    if not normalize:
        for protein in score_dict:
            score = round(float(score_dict[protein][0]))
            score_list.append(score)
        return score_list
    else:
        for protein in score_dict:
            score = float(score_dict[protein][0]) / float(score_dict[protein][1])
            score_list.append(score)
        return score_list
'''


def analise_membership(score_hash):
    group_hash = {}
    for protein in score_hash:
        protein = protein.strip(">")
        o_group = grap_output(protein, ortho_group_file)
        if o_group not in group_hash:
            group_hash[o_group] = 1
        else:
            group_hash[o_group] += 1
    return group_hash


def tp_vs_tn_plot(tp_list, tn_list, name):
    m_tp, std_tp = np.mean(tp_list), np.std(tp_list, ddof=1)
    if len(tn_list) >= 10:
        intersection = find_density_intersection(tp_list, tn_list)
        if intersection:
            single_inter = intersection
        else:
            single_inter = m_tp/2
            median = np.median(tp_list)/2
    else:
        single_inter = m_tp/2
        median = np.median(tp_list)/2
    intersection = find_density_intersection(tp_list, tn_list)
    sns.distplot(tp_list, hist=False, rug=True, color="r", label="PHMM intracluster scores")
    sns.distplot(tn_list, hist=False, rug=True, color="b", label="Next best BLAST PHMM scores")
    #plt.axvline(x=single_inter, linewidth=1, color="k", linestyle="dotted", label="mean")
    #plt.axvline(x=median, linewidth=1, color="k", linestyle="dashed", label="median")
    #if intersection:
    #    plt.axvline(x=intersection, linewidth=1, color="k", linestyle="solid", label="intersection point")
    plt.title("PHMM-Score Distribution: {}".format(name), size=16, weight="bold")
    plt.xlabel("PHMM scores", size=14)
    plt.ylabel("Density", size=14)
    plt.savefig(out_name)
    plt.close()


def plot_tn_group_belonging(true_negative_hash, cluster_name):
    tn_memberships = analise_membership(true_negative_hash)
    name_list = []
    count_list = []
    for name_x in tn_memberships:
        name_list.append(name_x)
        count_list.append(tn_memberships[name_x])
    ax = plt.subplot(111)
    width = 0.3
    bins = [x - (width/2) for x in range(1, len(name_list)+1)]
    ax.bar(bins, count_list, width=width)
    ax.set_xticks([x + (width/2) for x in range(1, len(name_list)+1)])
    ax.set_xticklabels(name_list, rotation=30, size=10, rotation_mode="anchor", ha="right", weight="bold")
    plt.title("Orthologous Groups of next best BLAST proteins of {}".format(cluster_name), size=12, weight="bold")
    plt.ylabel("Group Size", size=14)
    plt.savefig(out_name)
    plt.close()
    return tn_memberships


def walk_and_hash_input(input_dir):
    tn_pathes = {}
    tp_pathes = {}
    for subdir, dirs, files in os.walk(input_dir):
        for single_file in files:
            file_path = os.path.join(os.path.abspath(subdir), single_file)
            single_file = ".".join(single_file.strip().split(".")[0:-1]).split("_")
            file_belonging = single_file[-1]
            file_name = "_".join(single_file[:-1])
            #file_belonging = file_path.split("_")[-1]
            #file_name = single_file[-1]
            if file_belonging == "scores":
                tp_pathes[file_name] = file_path
            elif file_belonging == "TrueNegativeScores":
                tn_pathes[file_name] = file_path
            else:
                print("path Error")
                pass
    return tn_pathes, tp_pathes


if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    input_dir = args['--input_dir']
    output_dir = get_outdir(args['--out_dir'])
    ortho_group_file = args['--membership']

    tn_files, tp_files = walk_and_hash_input(input_dir)

    summary_membership = {}
    for name in tn_files:
        tp_file = tp_files[name]
        tn_file = tn_files[name]
        tp_hash, tp_scores = hash_score_files(tp_file)
        tn_hash, tn_scores = hash_score_files(tn_file)
        try:
            name = name.split(".")[-1]
            if "_" in name:
                name = name.split("_")[1]
        except Exception:
            pass
        out_name = os.path.join(output_dir, name + "_PHMMscores.pdf")
        if len(tn_scores) > 1 and len(tp_scores) > 1:
            tp_vs_tn_plot(tp_scores, tn_scores, name)
            if args["--membership"]:
                out_name = os.path.join(output_dir, name + "_membership.pdf")
                summary_membership[name] = plot_tn_group_belonging(tn_hash, name)
            print(name)
        else:
            print("{} not enough values to plot".format(name))
    if args["--membership"]:
        with open(os.path.join(output_dir,"NBB_membership_table.txt"), "w") as sum_f:
            sum_f.write("### Summary of Orthofinder groups which had next best BLAST hits to another cluster\n### first column = group name\tsecond column = count\n\n")
            for cluster_name, nbb_dict in summary_membership.items():
                sum_f.write("#{}\n".format(cluster_name))
                for group_name, count in nbb_dict.items():
                    sum_f.write("{}\t{}\n".format(group_name, str(count)))
                sum_f.write("\n")
