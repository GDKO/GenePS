#!/usr/bin/env python3
# Script takes a RAxmL file and a corresponding fasta file as input and checks
# whether all species show a one-to-one-relationship in the cluster.
# The Fasta header have to be in the format in style of >species_name.protein_name
# It's also possible to lower the perecent of species that need to be monophyletic: "-p" and to print a labeled tree: --tree
# if -o is specified, the longest isoform for an in-paralog, plus all singel species sequences will be written into a fasta file


"""
Usage: One_to_one_parser.py                        -r <FILE> -f <FILE> [-o <FILE>] [-p <INT>] [--tree] [--print_mono]

    Options:
        -h, --help                          show this screen.

        General
        -r, --raxml <FILE>                  raxml tree file (1-line)
        -f, --fasta <FILE>                  corresponding fasta file to the provided tree (header = species_name.protein_name)
        Optional
        -o, --output <FILE>                 output file name
        -p, --proportion <INT>              minimum percent of monophyletic nodes. (Default 100%)
        --tree                              print labeled ETE tree to the screen (label: species, sequences length)
        --print_mono                        if cluster is rejected, print all monophyletic nodes to the screen
"""

import sys
import os
from ete3 import Tree
from collections import Counter
from docopt import docopt


def hash_fasta(fasta_file):
    """takes path to fasta file and returns a dictionary in style of: >header:ADFASDF"""
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


def read_in_tree(tree_path):
    """reads_in_raxml tree file and returns tree in ETE format"""
    with open(tree_path) as t:
        for line in t:
            return Tree(line)


def create_label_dicts(fasta_hash):
    """Takes fasta_hash in style of header: seq. and Returns 2 dicts 1) protein_name : species; 2) protein_name : len(seq)
     !!!!! species name is inferred from fasta header!!!!
     The first part of the header should be the species name sperated by '.' and followed by the seq. name !!!"""
    protein_to_species = {}
    protein_to_length = {}
    for header, seq in fasta_hash.items():
        header = header.strip(">")
        if not len(header.split(".")) > 1:
            raise AttributeError("no species names found in fasta header! Wanted: >speciesName.proteinName")
        protein_to_species[header] = header.split(".")[0]
        protein_to_length[header] = len(seq)
    return protein_to_species, protein_to_length


def add_label_to_tree(tree, species_dict, length_dict):
    """permanently labels the ETE-tree with species-names in the first column and sequence-length in the second """
    for leaf in tree:
        leaf.add_features(species=species_dict.get(leaf.name, "none"))
        leaf.add_features(length=length_dict.get(leaf.name, "none"))


def split_species_into_paralogs_and_singles(protein_to_species_dict):
    in_paralog_species, single_species_proteins = [], []
    species_count = Counter(protein_to_species_dict.values())
    for protein, species in protein_to_species_dict.items():
        if species_count[species] > 1:
            in_paralog_species.append(species)
        else:
            single_species_proteins.append(protein)
    return set(in_paralog_species), set(single_species_proteins)


def check_monophyly(tree, set_species_label):
    """checks for monophyly for each label(species) in the label set and returns positive tree nodes stored in a dictionary"""
    mono_dict = {}
    for species in set_species_label:
        if tree.check_monophyly(values=[species], target_attr="species")[0] is True:
            mono_dict[species] = tree.search_nodes(species=species)
    return mono_dict
    # for node in t.get_monophyletic(values=[color], target_attr="color"):
    # print((node.get_ascii(attributes=["color", "name"], show_internal=False)))


def get_longest_in_paralog(tree, mono_dict, length_dict):
    """for each monophyletic groups it picks the longest sequence and returns hash in style of: species : protein_name"""
    species_to_longest_protein_name = {}
    for species, obj_list in mono_dict.items():
        for node in obj_list:
            if species not in species_to_longest_protein_name:
                species_to_longest_protein_name[species] = node.name
            elif length_dict[species_to_longest_protein_name[species]] < list(tree.get_cached_content(store_attr="length")[node])[0]:
                species_to_longest_protein_name[species] = node.name
            else:
                pass
    return set(species_to_longest_protein_name.values())


def write_hash_to_fasta(file_path, protein_header_set, fasta_hash):
    """writes fasta file and takes for every protein header in results hash the corresponding sequence from fasta_hash"""
    if not protein_header_set:
        return None
    with open(file_path, "w") as protein_f:
        for protein_header in protein_header_set:
            protein_header = ">{}".format(protein_header)
            protein_f.write("{}\n{}\n".format(protein_header, fasta_hash[protein_header]))
    return file_path


def subset_fasta_dict(old_fasta_dict, ete3_tree):
    # makes fasta file that only contains entries which are also present in the raxml tree
    new_dict = {}
    for node in ete_tree:
        header = ">" + node.name
        new_dict[header] = old_fasta_dict[header]
    return new_dict


if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    tree_path = args['--raxml']
    fasta_file = args['--fasta']
    outpath = args['--output']
    print_tree = args['--tree']
    print_mono = args["--print_mono"]
    proportion_monophyl = 100
    if args['--proportion']:
        proportion_monophyl = float(args['--proportion'])

    fasta_dict_input = hash_fasta(fasta_file)
    ete_tree = read_in_tree(tree_path)
    fasta_dict = subset_fasta_dict(fasta_dict_input, ete_tree)
    species_label, length_label = create_label_dicts(fasta_dict)
    paralog_species_set, single_species_protein_header = split_species_into_paralogs_and_singles(species_label)
    add_label_to_tree(ete_tree, species_label, length_label)
    if print_tree:
        print("\n# Tree = header, species label, sequence length")
        print(ete_tree.get_ascii(attributes=["name", "species", "length"], show_internal=True))
    if not paralog_species_set:
        print("\n[!] all proteins are from different species -> files remained unchanged \n")
        sys.exit()
    try:
        monophyl_dict = check_monophyly(ete_tree, paralog_species_set)
        # number monophyl species / number total species
        percent_mono = ((len(monophyl_dict) + len(single_species_protein_header)) /
                        (len(paralog_species_set) + len(single_species_protein_header))) * 100
    except ValueError:
        print("\n[!] RAxML tree does not match with the Fasta File! Please check if all entries in the protein file can be found in the RAxML file...(or maybe wrong file pair?)\n")
        sys.exit()
    if percent_mono == 0:
        print("\n[!] No 1-to-1-relationship\n")
        sys.exit()
    elif percent_mono >= proportion_monophyl:
        print("\n[+] OK: {} - {}% are monophyletic".format(os.path.split(tree_path)[-1], round(percent_mono), 2))
        if outpath:
            species_to_1_inparalog = get_longest_in_paralog(ete_tree, monophyl_dict, length_label)
            inparalog_plus_singles_header = species_to_1_inparalog.union(single_species_protein_header)
            print("\n[+] Results written to: ", write_hash_to_fasta(outpath, inparalog_plus_singles_header, fasta_dict))
    else:
        print("\n[!] {} rejected: only {}% are monophyletic\n".format(os.path.split(tree_path)[-1], round(percent_mono), 2))
        if print_mono:
            print("[+] monophyletic nodes :")
            for species, nodes in monophyl_dict.items():
                print(species, ":\t", ",".join([x.name for x in nodes]))
            for singles_header in single_species_protein_header:
                print (singles_header)
    print("\n")
