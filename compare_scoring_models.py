#!/usr/bin/env python3

"""
Usage: build_models.py                         -i <FILE> -r <INT> -p <INT> -o <STR> -l <INT> [-q <FILE>]

    Options:
        -h, --help                            show this screen.

        General
        -i, --input <DIR>                     protein fasta-file from one cluster
        -r, --iterations <DIR>                number of iterations to increase the data set size
        -p, --number_of_proteins <FILE>       number of proteins added to the data set per iteration
        -l, --repeat_limit <INT>              how often to repeat the entire computation to plot min, max curves
        -o, --output_path <STR>               name and path of the results graph in pdf format

"""

import matplotlib.pyplot as plt
import seaborn as sns
import tempfile
from docopt import docopt
from random import sample
from shared_code_box import run_cmd, tempdir, hash_fasta, generate_hmm, write_to_tempfile
from build_models import ScoreObject, generate_msa
sns.set(color_codes=True)


def add_random_proteins(number):
    count = 0
    while count < number:
        random_header = sample(cluster_hash.keys(), 1)[0]
        if random_header not in seq_header_list:
            seq_header_list.append(random_header)
            count += 1


def iterative_scoring_model():
    with tempdir() as tmp:
        intra_score = ScoreObject(cluster_hash, seq_header_list, "cele_scoring", tmp)
        score_list = intra_score.iterative_score_computation()
    return score_list


def many_hmm_scores(hmm_file, query_file, test_this_function=None):
    if test_this_function is None:
        command = ["hmmsearch", "--noali", hmm_file, query_file]
    else:
        command = test_this_function
    read_count = 0
    score_list = []
    for line in run_cmd(command=command, wait=False):
        if "E-value" in line or read_count == 1:
            read_count += 1
        elif read_count == 2:
            line = line.strip("\n").split()     # line[8] is protein name
            if len(line) > 0:
                try:
                    score_list.append(float(line[1]))
                except ValueError:
                    print("VALUE ERROR", line)
                    return score_list
            else:
                return score_list


def bulk_hmm_scoring_model():
    with tempfile.NamedTemporaryFile() as tmp:
        seq_list = []
        for header in seq_header_list:
            seq_list.append(header)
            seq_list.append(cluster_hash[header])
        write_to_tempfile(tmp.name, "\n".join(seq_list))
        msa = "\n".join(generate_msa(tmp.name))
        with tempfile.NamedTemporaryFile() as msa_tmp:
            write_to_tempfile(msa_tmp.name, msa)
            with tempfile.NamedTemporaryFile() as tmp_hmm:
                hmm = generate_hmm(tmp_hmm.name, msa_tmp.name)
                score_list = many_hmm_scores(hmm, tmp.name)
                return score_list


def compare_models(runs, interval):
    average_list_single_hmm = []
    average_list_multiple_hmm = []
    counter = 1
    while counter <= runs:
        add_random_proteins(interval)
        try:
            print("\nBulk scoring")
            single_hmm_scores = bulk_hmm_scoring_model()
            single_hmm_mean = round(sum(single_hmm_scores) / len(single_hmm_scores))
        except TypeError:
            single_hmm_mean = None
        try:
            print("\nIterative scoring\n")
            multiple_hmm_scores = iterative_scoring_model()
            multiple_hmm_mean = round(sum(multiple_hmm_scores) / len(multiple_hmm_scores))
        except TypeError:
            multiple_hmm_mean = None
        average_list_multiple_hmm.append(multiple_hmm_mean)
        average_list_single_hmm.append(single_hmm_mean)
        counter += 1
    return average_list_single_hmm, average_list_multiple_hmm


def make_array(target_list):
    sum_list = []
    for idx in range(0, len(target_list[0])):
        pos_list = []
        for single_list in target_list:
            pos_list.append(single_list[idx])
        sum_list.append(pos_list)
    return sum_list


def mean_list(x):
    return round(sum(x) / len(x))


def apply_function_on_list_in_array(target_list, function):
    results = []
    for sub_list in target_list:
        try:
            results.append(function(sub_list))
        except TypeError:
            results.append(0)
    return results


if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    cluster_file = args['--input']
    number_prot = int(args['--number_of_proteins'])
    iterations = int(args['--iterations'])
    repeat_limit = int(args['--repeat_limit'])
    out_dir = args['--output_path']
    out_file = open(out_dir + ".txt", "w")
    # iterations = how often gets the dataset enlarged
    # number_prot = number of proteins added to the dataset per iteration
    # iterations * numbers = total amount of protein seq
    # per iteration get average scores for each model till iter * number_prot is reached
    # results in an "average-list" for both models

    # for each repeat step, repeat this computation form above
    # for each repeat step one average_list
    # average the average_lists
    # plot average of average and min of average and max of average for both models

    cluster_hash = hash_fasta(cluster_file)
    x_axis = [number_prot*x for x in range(1, iterations+1)]

    repeat_step = 1
    seq_header_list = []
    bulk_av_list = []
    iter_av_list = []
    while repeat_step <= repeat_limit:
        av_single_hmm, av_multiple_hmm = compare_models(iterations, number_prot)
        print("{} round of {}". format(str(repeat_step), str(repeat_limit)))
        print(av_single_hmm)
        print(av_multiple_hmm)
        bulk_av_list.append(av_single_hmm)
        iter_av_list.append(av_multiple_hmm)
        seq_header_list = []
        repeat_step += 1

    bulk_results = make_array(bulk_av_list)
    [out_file.write(str(x) + "\n") for x in bulk_results]
    bulk_min = apply_function_on_list_in_array(bulk_results, min)
    bulk_max = apply_function_on_list_in_array(bulk_results, max)
    bulk_mean = apply_function_on_list_in_array(bulk_results, mean_list)
    out_file.write("######iterational results\n")
    iter_results = make_array(iter_av_list)
    iter_min = apply_function_on_list_in_array(iter_results, min)
    iter_max = apply_function_on_list_in_array(iter_results, max)
    iter_mean = apply_function_on_list_in_array(iter_results, mean_list)
    [out_file.write(str(x) + "\n") for x in iter_results]
    out_file.close()

    plt.plot(x_axis, bulk_mean, color="b", label="Bulk HMM-Scoring")
    plt.fill_between(x_axis, bulk_min, bulk_max, facecolor="LightSkyBlue", alpha=0.5)
    plt.plot(x_axis, iter_mean, color="r", label="Iterative HMM-Scoring")
    plt.fill_between(x_axis, iter_min, iter_max, facecolor="Salmon", alpha=0.5)

    plt.title("Comparison of Scoring Models", size=18, weight="bold")
    plt.ylabel("Score Size", size=14)
    plt.xlabel("Number of Protein Sequences", size=14)
    plt.legend()
    plt.savefig(out_dir + ".pdf")
    plt.show()

