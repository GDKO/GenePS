# GenePS
### Gene Positioning System

UNDER CONSTRUCTION  !! :)

GenePS is a protein cluster based, local gene finder and aims to find genes in newly sequenced genomes in the absence
of any expression data from the species itself. GenePS merely requires a collection protein sequences from related species as input files.

SHORT MANUAL

Download all files of GenePS from Github. Most of the files are libraries but need to be present at the current stage
of development. Furthermore the following programs and python packages are required to run GenePS:

Programs:
- HMMER
- MAFFT
- BLAST
- Exonerate
- Trimal

Python package:
- (Python 3)
- docopt
- matplotlib.pyplot
- seaborn
- numpy
- scipy
- sklearn

GenePS consists of two separate scripts: "build_models.py" and "use_models.py". The "build_models.py" script generates
an output directory with profile HMMs and parameter files which can then, at any time, be deployed by "use_models.py"
to carry out the actual gene finding.


1) make_GenePS.py
~~~
Usage: build_models.py                         -i <DIR> -o <DIR> [-f <FILE>] [--keep] [--subset] [--einsi]

    Options:
        -h, --help                            show this screen.

        General
        -i, --input <DIR>                     either single fasta file or folder with only fasta files
        -o, --output <DIR>                    creates a directory for all output files
        -f, --orthofinder_files <FILE>        4 lines file in style of: blast_dir=DIR speciesIDs=file.txt proteins=protein.fa sequenceIDs=file.txt
        --keep                                command to safe intermediate files
        --subset                              clusters all filtered proteins by length and outputs a length representative subset instead of all filtered proteins
        --einsi                               changes the default simple progressive method of mafft to the E-INS-i algorithm
~~~
The input can be either a single fasta file containing raw protein sequences from one gene family or
a whole directory with many protein files. It is also possible to group the protein files in folders. Fasta files in a folder
are then labeled by the folder's name they are kept in during the build models process. For instance: names of fasta files could be gene names and names
of folders could be species names.
The output will be a profile HMM, a next best BLAST profile HMM and a filtered protein file for every input cluster.
Furthermore and additional fasta file with consensus sequences as well as a parameter.GenePS file will be generate, one per folder (storing information for all files in this folder)


2) use_models.py
~~~
Usage: use_models.py                          -i <DIR> -g <FILE> [-c <INT>] [-o <DIR>] [--keep] [--verbose] [--frag] [--quick]

    Options:
        -h, --help                            show this screen.

        General
        -i, --use_models_input <DIR>          Directory with results from build_models.py
        -g, --genome <FILE|LIST>              Single genome fasta-file or a many lines ".txt"-file in style of: name=path/to/genome
        Optional
        -c, --coverage_filer <INT>            Minimal aligned length of a Blast query to the target genome in % (used to filter Blast hits)
        -o, --out_dir <DIR>                   Directory for the output files and folders (default: same as input directory)
        --quick                               Omits exhaustive prediction refinements through and additional exonerate -E yes run (faster but less accurate)
        --frag                                If enabled, a length filter will be applied to remove potentially fragmented predictions
        --keep                                Keeps intermediate files (Blast output, merged regions, exonerate output)
        --verbose                             Prints progress details to the screen


~~~

"use_models.py" takes the previous output directory as input. The second demanded input is the genome in fasta format. Alternatively,
multiple genomes can be specified by creating a text if in style of:

genome_name=path/to/genome

...

genome_name_N=path/to/genome_N

Please note that there are no spaces between the equal sign and the path variable.
This text file can then be given to GenePS by -g genomes.txt.

optional arguments:
coverage_filter: GenePS uses TBLASTN align the consensus sequences against the genome define candidate region. Regions
are only considered for further downstream analysis if a certain percent of the query length aligns to this region. The
default minimum for query coverage is 30%. A higher threshold could exclude false positive regions and save run time.
If the threshold is too low a true positive region might be missed.

out_dir: specify an output directory. If nothing is specified a "Predictions"-folder will be created
in the same directory as the "build_modles.py"-output.

quick: GenePS first aligns all protein sequences against a candidate region to find the protein sequence within a cluster
that is most likely to give the best prediction for this particular region. Alignments and predictions in GenePS are
generate by exonerate. As exonerate is a time consuming program, this first "query selection" step is carried out in the "p2g -E no" mode.
In the default mode of "use_models.py" a second refinement step follows. This is done by re-aligning the previously identified,
most promising sequence to the same region again, but now in the "p2g -E yes" mode (exhaustive). In order save run time the second
refinement step can be omitted with the --quick flag. It will speed up the run dramatically but often on the costs of accuracy.

frag: Activates a fragmentation filter. Default off. If on, every protein predicted by exonerate will have to be within
a certain sequence length range in order to pass the filter of GenePS. Predictions outside this length range will be
considered as fragmented predictions. The lenght range was computed during model building and is stored in the .GenePS file.

--keep: keeps intermediate results file. Blast files, merged HSP files and exonerate intermediate files.

--verbose: Prints progress statements to the screen (e.g. about passed predictions or even filtered prediction).


The final output of "use_models" will be gff files, protein files and cds files for predictions GenePS considers to be "valid"
but also for all filtered predictions. Furthermore a summary report will be written to the output directory as well as a LOG file.
