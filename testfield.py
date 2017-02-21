#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.signal import argrelmax
from sklearn.neighbors.kde import KernelDensity
from sklearn import cluster
import statsmodels.sandbox
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy import stats
from plot_PHMM_scores import walk_and_hash_input, hash_score_files
from build_models import find_density_intersection
import random
import seaborn as sns
import pandas
from tkinter import *

dir_default_scores = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/mafft_default_lethal_clusters/intermediate_files"
dir_einsi = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/mafft_einsi_lethal_clusters/intermediate_files"
tn_files_default, tp_files_default = walk_and_hash_input(dir_default_scores)
tn_files_einsi, tp_files_einsi = walk_and_hash_input(dir_einsi)
pval_list = []
means_default = []
means_einsi = []
means_array =[]
std_default, std_einsi, intersect_diff= [], [], []
for name in tn_files_default:
    # mafft default parameters and variables
    tp_file_default = tp_files_default[name]
    tn_file_default = tn_files_default[name]
    tp_hash_default, tp_scores_default = hash_score_files(tp_file_default)
    tn_hash_default, tn_scores_default = hash_score_files(tn_file_default)
    means_default.append(np.mean(tp_scores_default))
    std_default.append(np.std(tp_scores_default, ddof=1))
    # mafft einsi paramters and variables
    tp_file_einsi = tp_files_einsi[name]
    tn_file_einsi = tn_files_einsi[name]
    tp_hash_einsi, tp_scores_einsi = hash_score_files(tp_file_einsi)
    tn_hash_einsi, tn_scores_einsi = hash_score_files(tn_file_einsi)
    means_einsi.append(np.mean(tp_scores_einsi))
    std_default.append(np.std(tp_scores_einsi, ddof=1))
    # diff in intersection
    if len(tp_scores_einsi) >= 10 and len(tn_scores_einsi) >= 10:
        intersect_default = find_density_intersection(tp_scores_default, tn_scores_default)
        intersect_einsi = find_density_intersection(tp_scores_einsi, tn_scores_einsi)
    else:
        intersect_default = np.median(tp_scores_default) / 2
        intersect_einsi = np.median(tp_scores_einsi) / 2
        print("median")
    print(abs(intersect_einsi-intersect_default), name)
    intersect_diff.append(str(round((abs(intersect_einsi-intersect_default)/intersect_default)*100, 1)) + "%")
    # append wilcoxon p values
    sample_siz = min([len(tp_scores_default), len(tp_scores_einsi)])
    wilcoxon = stats.wilcoxon(random.sample(tp_scores_default, sample_siz), random.sample(tp_scores_einsi, sample_siz))
    pval_list.append(wilcoxon[1])

print(means_default)
print(means_einsi)

# make multiple testing correction
print(multipletests(pval_list, method="s"))



# plot
sns.set_style("darkgrid")
plt.axis([0, len(means_default)+1, min(means_default)-60, max(means_einsi)+50])

line1 = plt.Line2D(range(1), range(1), color="whitesmoke", marker='o', markerfacecolor="red")
line2 = plt.Line2D(range(1), range(1), color="whitesmoke", marker='s', markerfacecolor="blue")

legend_properties = {"weight":"bold", "size":12}
plt.legend((line1, line2), ('MAFFT: default','MAFFT: einsi'), loc='upper right', prop=legend_properties)

plt.plot(range(1,len(means_default)+1), means_einsi, "bs", color="blue")
plt.plot(range(1,len(means_default)+1), means_default, "ro", color="red")

count = 0 # add percent of intersect difference
for label in intersect_diff:
    y = means_default[count] - 35
    count += 1
    plt.text(count-0.2, y, label, color="black", size=12, weight="bold")

plt.ylabel('PHMM scores', size=16, weight="bold")
plt.xlabel('Cluster', size=16, weight="bold")
plt.title("Effect of MAFFT modes on pHMM scores", size=16, weight="bold")
plt.show()
