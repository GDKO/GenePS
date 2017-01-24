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


dir_default_scores = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/mafft_option_test/mafft_default_test/intermediate_files"
dir_einsi = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/mafft_option_test/einsi_test/intermediate_files"
tn_files_default, tp_files_default = walk_and_hash_input(dir_default_scores)
tn_files_einsi, tp_files_einsi = walk_and_hash_input(dir_einsi)
pval_list = []
means_default = []
means_einsi = []
means_array =[]
std_default, std_einsi = [], []
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
    intersect_default = find_density_intersection(tp_scores_default, tn_scores_default)
    intersect_einsi = find_density_intersection(tp_scores_einsi, tn_scores_einsi)
    print(abs(intersect_einsi-intersect_default), name)
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

plt.axis([-1,6, 400,1600])
p1 = plt.Rectangle((0, 0), 0.1, 0.1, fc="crimson")
p2 = plt.Rectangle((0, 0), 0.1, 0.1, fc="burlywood")

plt.legend((p1, p2), ('category 1','category 2'), loc='upper left')

plt.plot(means_einsi, "bs", color="burlywood")
plt.plot(means_default, "ro", color="crimson")

plt.ylabel('This is your y-axis')
plt.xlabel('This is my x-axis')
plt.title("This is our title")
plt.show()
