#!/usr/bin/env python3
import subprocess
import shutil
import sys
import contextlib
import os
import tempfile as tmp
from collections import defaultdict


def hash_fasta(fasta_file):
    """takes sting to fasta file and returns a dictionary in style of: >header:[ADFASDF]"""
    try:
        fasta = {}
        active_sequence_name = ""
        with open(fasta_file) as file_one:
            for line in file_one:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    active_sequence_name = line.split(" ")[0]
                    if active_sequence_name not in fasta:
                        fasta[active_sequence_name] = []
                    continue
                sequence = line
                fasta[active_sequence_name].append(sequence)
        file_one.close()
        return {header: "".join(seq) for header, seq in fasta.items()}
    except KeyError:
        return None


def write_hash_to_fasta(file_path, dictionary, line_style="{}\n{}\n"):
    """writes a hash to a file in the given line_style. Line style default writes key and value to separate lines.
    If value is a list it joins elements without spaces."""
    if not dictionary:
        return None
    try:
        with open(file_path, "w") as hash_file:
            for header, value in dictionary.items():
                hash_file.write(line_style.format(header, value))
        return file_path
    except IndexError or KeyError or AttributeError or TypeError:
        return None


def get_outdir(out_directory, add_dir=""):
    """generates output directory in case it does not exist."""
    if type(out_directory) != str:
        print("\t[!] {} is NOT a directory! Please specify an output directory\n".format(out_directory))
        sys.exit()
    elif os.path.isfile(out_directory):
        print("\t[!] {} is a File! Please specify an output directory\n".format(out_directory))
        sys.exit()
    elif not os.path.exists(os.path.join(out_directory, add_dir)):
        os.mkdir(os.path.join(out_directory, add_dir))
        return os.path.abspath(os.path.join(out_directory, add_dir))
    else:
        return os.path.abspath(os.path.join(out_directory, add_dir))


def generate_hmm(hmm_path, msa_path):
    """builds hidden markov model from MSA"""
    command = ["hmmbuild", hmm_path, msa_path]
    run_cmd(command=command, wait=True)
    return hmm_path


def get_phmm_score(hmm_file, query_file):
    """aligns query-fasta-file against HMM and returns a score hash in style of >header:score.
    The score is defined by: sum(score domains)/Ndomains * sum(DOMAINend-DOMAINstart/query_lengt, ...)"""
    with tmp.NamedTemporaryFile() as domtblout:
        command = ["hmmsearch", "--domtblout", domtblout.name, "--noali", hmm_file, query_file]
        run_cmd(command=command, wait=True)
        final_score_hash = parse_hmmer_domain_table(domtblout.name)
    return final_score_hash


def parse_hmmer_domain_table(hmmer_table):
    score_dict = {}
    coverage_dict = defaultdict(list)
    for line in open(hmmer_table):
        if not line.startswith("#"):
            line = line.strip("\n").split()
            if ">" + line[0] not in coverage_dict:
                score_dict[">" + line[0]] = float(line[7])/float(line[2])
                coverage_dict[">" + line[0]].append([float(line[17]), float(line[18])])
            elif coverage_dict[">" + line[0]][-1][1] >= float(line[17]):
                coverage_dict[">" + line[0]][-1][1] = float(line[18])
            else:
                coverage_dict[">" + line[0]].append([float(line[17]), float(line[18])])
    return {header: round(score_dict[header] * (sum([x[1]-x[0] for x in coverage]))) for header, coverage in coverage_dict.items()}


def get_consensus(hmm_file):
    """infers a header-less consensus sequence from HMM"""
    command = "hmmemit -c " + hmm_file
    cons_list = []
    for line in run_cmd(command=command, wait=False):
        if not line.startswith(">"):
            cons_list.append(line.strip("\n"))
    return "".join(cons_list)


def print_progress(iteration, total, prefix='', suffix='', decimals=1, bar_length=100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    format_str = "{0:." + str(decimals) + "f}"
    percents = format_str.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = '█' * filled_length + '-' * (bar_length - filled_length)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.write('\r')
    if iteration == total:
        sys.stdout.write("\n")
    sys.stdout.flush()


def write_to_tempfile(tmp_name, string):
    """takes tmp_name, opens it as a real file and writes string into it"""
    new_file = open(tmp_name, "w")
    new_file.write(string)
    new_file.seek(0)
    new_file.close()


'''
def find_intersection_normal_distribution(m_p, m_n, std_p, std_n):
    """find the intersection between to distributions assuming they are normally distributed. If > 1 insersection points
    it returns the biggest inter.-point which is smaller then the mean of the true positive distribution"""
    a = 1/(2 * std_p ** 2) - 1 / (2 * std_n ** 2)
    b = m_n / (std_n ** 2) - m_p / (std_p ** 2)
    c = m_p ** 2 / (2 * std_p ** 2) - m_n ** 2 / (2 * std_n ** 2) - log(std_n / (std_p * 2))
    intersection = roots([a, b, c])
    sorted_no_negative = sorted([round(x, 3) for x in intersection if not x < 0], reverse=True)
    for point in sorted_no_negative:
        if point < m_p + std_p:
            return point
    return None
'''


def run_cmd(**kwargs):
    cmd = kwargs['command']
    if type(cmd) != list:
        cmd = kwargs['command'].split() # sanitation
    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,
                         universal_newlines=True,
                         bufsize=-1)
    wait = kwargs.get('wait', False)
    if wait:
        stream_data = p.communicate()[0]
        p.wait()
        if p.returncode == 0:
            pass
    else:
        return iter(p.stdout.readline, '')


def which(program):
    if shutil.which(program):
        return program
    else:
        return False


def check_programs(*arg):
    error_list = []
    for program in arg:
        if which(program) is False:
            error_list.append("\t[!] {} not found! Please install and add to it PATH variable".format(program))
    if error_list:
        print("\n".join(error_list))


@contextlib.contextmanager
def cd(newdir, cleanup=lambda: True):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)
        cleanup()


@contextlib.contextmanager
def tempdir():
    dirpath = tmp.mkdtemp()

    def cleanup():
        shutil.rmtree(dirpath)
    with cd(dirpath, cleanup):
        yield dirpath


