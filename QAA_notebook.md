# 2024 Bi623 RNA-seq Quality Assessment Assignment Lab Notebook

 ------------------
| Sunday 240825    |
 ------------------

## Part 1 – Read quality score distributions

"The objectives of this assignment are to use existing tools for quality assessment and adaptor trimming, compare the quality assessments to those from your own software, and to demonstrate your ability to summarize other important information about this RNA-Seq data set in a high-level report. That is, you should create a cohesive, well written report for your "PI" about what you've learned about/from your data."

We've each been assigned two datasets (each containing R1 and R2) to work on during this project. We are not to move, copy, or unzip these data!

My pairs: Amelia  2_2B_control_S2_L008    19_3F_fox_S14_L008

Their paths are:
```
/projects/bgmp/shared/2017_sequencing/demultiplexed/2_2B_control_S2_L008_R1_001.fastq.gz
/projects/bgmp/shared/2017_sequencing/demultiplexed/2_2B_control_S2_L008_R2_001.fastq.gz

/projects/bgmp/shared/2017_sequencing/demultiplexed/19_3F_fox_S14_L008_R1_001.fastq.gz
/projects/bgmp/shared/2017_sequencing/demultiplexed/19_3F_fox_S14_L008_R2_001.fastq.gz
```

Creating an environment called QAA for this project:

```
$ conda create -n QAA python=3.12 FastQC
(QAA) [02:32:18 /projects/bgmp/shared/2017_sequencing/demultiplexed] adayton:n0353.talapas.uoregon.edu
$ fastqc --version
FastQC v0.12.1
```

We're good - that's the version I want. Time to figure out how to use fastqc on the CLI so I can do this assignment.

I need to: "produce plots of the per-base quality score distributions for R1 and R2 reads. Also, produce plots of the per-base N content, and comment on whether or not they are consistent with the quality score plots."

Here's a helpful web page: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

And the help option:

```
$ fastqc --help

            FastQC - A high throughput sequence QC analysis tool

SYNOPSIS

        fastqc seqfile1 seqfile2 .. seqfileN

    fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] 
           [-c contaminant file] seqfile1 .. seqfileN

DESCRIPTION

    FastQC reads a set of sequence files and produces from each one a quality
    control report consisting of a number of different modules, each one of 
    which will help to identify a different potential type of problem in your
    data.
    
    If no files to process are specified on the command line then the program
    will start as an interactive graphical application.  If files are provided
    on the command line then the program will run with no user interaction
    required.  In this mode it is suitable for inclusion into a standardised
    analysis pipeline.
    
    The options for the program as as follows:
    
    -h --help       Print this help file and exit
    
    -v --version    Print the version of the program and exit
    
    -o --outdir     Create all output files in the specified output directory.
                    Please note that this directory must exist as the program
                    will not create it.  If this option is not set then the 
                    output file for each sequence file is created in the same
                    directory as the sequence file which was processed.
                    
    --casava        Files come from raw casava output. Files in the same sample
                    group (differing only by the group number) will be analysed
                    as a set rather than individually. Sequences with the filter
                    flag set in the header will be excluded from the analysis.
                    Files must have the same names given to them by casava
                    (including being gzipped and ending with .gz) otherwise they
                    won't be grouped together correctly.
                    
    --nano          Files come from nanopore sequences and are in fast5 format. In
                    this mode you can pass in directories to process and the program
                    will take in all fast5 files within those directories and produce
                    a single output file from the sequences found in all files.                    
                    
    --nofilter      If running with --casava then don't remove read flagged by
                    casava as poor quality when performing the QC analysis.
                   
    --extract       If set then the zipped output file will be uncompressed in
                    the same directory after it has been created. If --delete is 
                    also specified then the zip file will be removed after the 
                    contents are unzipped. 
                    
    -j --java       Provides the full path to the java binary you want to use to
                    launch fastqc. If not supplied then java is assumed to be in
                    your path.
                   
    --noextract     Do not uncompress the output file after creating it.  You
                    should set this option if you do not wish to uncompress
                    the output when running in non-interactive mode.
                    
    --nogroup       Disable grouping of bases for reads >50bp. All reports will
                    show data for every base in the read.  WARNING: Using this
                    option will cause fastqc to crash and burn if you use it on
                    really long reads, and your plots may end up a ridiculous size.
                    You have been warned!
                    
    --min_length    Sets an artificial lower limit on the length of the sequence
                    to be shown in the report.  As long as you set this to a value
                    greater or equal to your longest read length then this will be
                    the sequence length used to create your read groups.  This can
                    be useful for making directly comaparable statistics from 
                    datasets with somewhat variable read lengths.

    --dup_length    Sets a length to which the sequences will be truncated when 
                    defining them to be duplicates, affecting the duplication and
                    overrepresented sequences plot.  This can be useful if you have
                    long reads with higher levels of miscalls, or contamination with
                    adapter dimers containing UMI sequences.

                    
    -f --format     Bypasses the normal sequence file format detection and
                    forces the program to use the specified format.  Valid
                    formats are bam,sam,bam_mapped,sam_mapped and fastq
                    

    --memory        Sets the base amount of memory, in Megabytes, used to process 
                    each file.  Defaults to 512MB.  You may need to increase this if
                    you have a file with very long sequences in it.
                
    --svg           Save the graphs in the report in SVG format.

    -t --threads    Specifies the number of files which can be processed
                    simultaneously.  Each thread will be allocated 250MB of
                    memory so you shouldn't run more threads than your
                    available memory will cope with, and not more than
                    6 threads on a 32 bit machine
                  
    -c              Specifies a non-default file which contains the list of
    --contaminants  contaminants to screen overrepresented sequences against.
                    The file must contain sets of named contaminants in the
                    form name[tab]sequence.  Lines prefixed with a hash will
                    be ignored.

    -a              Specifies a non-default file which contains the list of
    --adapters      adapter sequences which will be explicity searched against
                    the library. The file must contain sets of named adapters
                    in the form name[tab]sequence.  Lines prefixed with a hash
                    will be ignored.
                    
    -l              Specifies a non-default file which contains a set of criteria
    --limits        which will be used to determine the warn/error limits for the
                    various modules.  This file can also be used to selectively 
                    remove some modules from the output all together.  The format
                    needs to mirror the default limits.txt file found in the
                    Configuration folder.
                    
   -k --kmers       Specifies the length of Kmer to look for in the Kmer content
                    module. Specified Kmer length must be between 2 and 10. Default
                    length is 7 if not specified.
                    
   -q --quiet       Suppress all progress messages on stdout and only report errors.
   
   -d --dir         Selects a directory to be used for temporary files written when
                    generating report images. Defaults to system temp directory if
                    not specified.
                    
BUGS

    Any bugs in fastqc should be reported either to simon.andrews@babraham.ac.uk
    or in www.bioinformatics.babraham.ac.uk/bugzilla/
```

Going to try just running this I guess? I want the -o option to go to my "fastqc_out" directory.

```
(QAA) [02:44:47 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0353.talapas.uoregon.edu
$ fastqc -o fastqc_out/ --extract --delete /projects/bgmp/shared/2017_sequencing/demultiplexed/2_2B_control_S2_L008_R1_001.fastq.gz /projects/bgmp/shar
ed/2017_sequencing/demultiplexed/2_2B_control_S2_L008_R2_001.fastq.gz
```

OKAY! this gave me some good looking plots in my fastqc_out. I renamed and moved the per base quality score and per base N content. I also deleted all the other info that I don't need at the moment because it kinda clogged stuff up. This was really quick - took about a minute; I can rerun if need be. I'm doing it for the other samples now.

```
(QAA) [02:47:19 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0353.talapas.uoregon.edu
$ fastqc -o fastqc_out/ --extract --delete /projects/bgmp/shared/2017_sequencing/demultiplexed/19_3F_fox_S14_L008_R1_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/19_3F_fox_S14_L008_R2_001.fastq.gz
```

Moved, renamed, deleted other files for these too.

Running my perbasehist.py script on the other files. This is gonna take a while as it usually takes around an hour for each plot to be generated.

When it's done I need to compare/contrast my plots and a buncha info about them. 

```
	Command being timed: "./perbasehist.py -f /projects/bgmp/shared/2017_sequencing/demultiplexed/2_2B_control_S2_L008_R1_001.fastq.gz -l 101 -u True -y 42 -o ./part1_plots/perbasehist.py_plots/2_2B_control_S2_L008_R1_001_perbasehist.png"
	User time (seconds): 81.35
	System time (seconds): 0.22
	Percent of CPU this job got: 93%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:26.87
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 64196
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 35541
	Voluntary context switches: 2590
	Involuntary context switches: 190
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

	Command being timed: "./perbasehist.py -f /projects/bgmp/shared/2017_sequencing/demultiplexed/2_2B_control_S2_L008_R2_001.fastq.gz -l 101 -u True -y 42 -o ./part1_plots/perbasehist.py_plots/2_2B_control_S2_L008_R2_001_perbasehist.png"
	User time (seconds): 81.77
	System time (seconds): 0.14
	Percent of CPU this job got: 100%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:21.25
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 65940
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 36948
	Voluntary context switches: 524
	Involuntary context switches: 65
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

	Command being timed: "./perbasehist.py -f /projects/bgmp/shared/2017_sequencing/demultiplexed/19_3F_fox_S14_L008_R1_001.fastq.gz -l 101 -u True -y 42 -o ./part1_plots/perbasehist.py_plots/19_3F_fox_S14_L008_R1_001_perbasehist.png"
	User time (seconds): 225.74
	System time (seconds): 0.26
	Percent of CPU this job got: 100%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 3:45.81
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 67876
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 41720
	Voluntary context switches: 636
	Involuntary context switches: 97
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
    
	Command being timed: "./perbasehist.py -f /projects/bgmp/shared/2017_sequencing/demultiplexed/19_3F_fox_S14_L008_R1_001.fastq.gz -l 101 -u True -y 42 -o ./part1_plots/perbasehist.py_plots/19_3F_fox_S14_L008_R2_001_perbasehist.png"
	User time (seconds): 226.05
	System time (seconds): 0.25
	Percent of CPU this job got: 100%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 3:45.95
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 62200
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 41948
	Voluntary context switches: 521
	Involuntary context switches: 98
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

 ------------------
| Tuesday 240827   |
 ------------------

Checked on the commands and stuck the slurm output above in Sunday's entry.

Overall, looks like each took a few minutes... I should rerun fastqc in a slurm script to save the data...

```
	Command being timed: "fastqc -o fastqc_out/ --extract --delete /projects/bgmp/shared/2017_sequencing/demultiplexed/2_2B_control_S2_L008_R1_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/2_2B_control_S2_L008_R2_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/19_3F_fox_S14_L008_R1_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/19_3F_fox_S14_L008_R2_001.fastq.gz"
	User time (seconds): 174.46
	System time (seconds): 7.84
	Percent of CPU this job got: 98%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 3:05.61
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 402412
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 122841
	Voluntary context switches: 34913
	Involuntary context switches: 293
	Swaps: 0
	File system inputs: 0
	File system outputs: 9576
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

Ok cool. Finished part 1. Starting part 2!

## Part 2 – Adaptor trimming comparison

Installing cutadapt and trimmomatic into my QAA environment.

```
mamba install cutadapt
mamba install trimmomatic
```

```
(QAA) [11:22:39 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0349.talapas.uoregon.edu
$ cutadapt --version
4.9
(QAA) [11:22:43 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0349.talapas.uoregon.edu
$ trimmomatic -version
0.39
```

Looks good! Those are the versions I need.

So first I need to use cutadapt to properly trim adapter sequences using the default settings. I need to report what proportion of reads (both R1 and R2) were trimmed. I also need to figure out what the adapters are!

Looking at the fastqc plots they appear mostly toward the end of the control sequences. I'm going to start there.

I did this to grab just some sequences to examine for similarities...

```
$ zcat /projects/bgmp/shared/2017_sequencing/demultiplexed/2_2B_control_S2_L008_R1_001.fastq.gz | grep -A 1 "^@" | grep -v "^@" | grep -v "^--$" | head
CNCTGTTTGATGACAGAGCAGGGTCTTGCTGAACTGCACCTACTGGAGTAGGCAACTGTTTACCAATCCAGTCATATTCATAATCAAACATATACCCTTTT
GNCTACTACATAGTATGTATCGTGAAGCACGATGTCAAGGGATGAGTTGGATAAAACAATTCCGGTTAGACCACCAACTGTAAATAAGAAAATAAAGCCTA
CNCACGTTGTCTGCTCCCACAATGAAGCATTTTGGATAATCATCCAAAAGTTGGATGATCTTGAGGAAGTAGTTGGACTTCCAGGTCGCCCTGTCTTCCCT
GNCCCGGAGGCTGCGACTGGAGAGACTTGTAAAGGCGGCTGGAAGCGAGTACAGTGGGAATCCCCGGGTCTGTTCCGTTCAAACACTGTTGAAGAGATCGG
TNGCTTTGCATTTGAAAGGTAAGGAGTTTCCAGGAAGAGAAGATCCACCTTCAATCCCAAGAGTCCAAAGAGGGTTCTCTGTGATACCCGGCTTGCAGAGA
TNGTGATATTCTCTGTTGTGACTTGTGGTGATTTGCAAAGGACTCTAGTTCTGTCTTCCGTTATAGGTTTATCTTTGGGGTTTGGTGTTTCAAAGTGTTCC
GNGAGCTTATGGATCTAGAGGAGTCTATGATGGAGAAATGAAGGCATGAGAATCTGTTCACATCTTGAAAGATCGGAAGAGCACACGTCTGAACTCCAGTC
CNCCGGATGCCAAAGGACCTGACAGTGTATCTAGCTTTGGAGAACACCGGGGTCTGGCCTGTGAGCTGCTCCAACACCTTGGCTGCCCGGGTCAGTCAGAT
ANGCCCCAAACCCAAACAACACACACACACACACACACACACACACACACACACACAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATCGATAT
CNCTAGGGTGGCCTGCTGCTCTTCGGCTCCCAGCCTCTGATAATTTGCCTTCAGTGTTTGACAGGCTGTTGCCTCCTATCTTTTGCTGCTTCTTGGGGCTA
(QAA) [11:26:45 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0349.talapas.uoregon.edu
$ zcat /projects/bgmp/shared/2017_sequencing/demultiplexed/2_2B_control_S2_L008_R2_001.fastq.gz | grep -A 1 "^@" | grep -v "^@" | grep -v "^--$" | head
NTCTAGGAAAAGAGCAGAGTAGAAGAGATGATTTAGAAGCCTTAGGCCACATGTTCATGTACTTCCTGAGAGGCAGTCTTCCTTGGCAAGGCTTAAAGGCT
NGGTGTCAAAGTATTTAGCTGACTTGCAACCCTACACGGAGGTAATATTAAATGATCTCCAGCTATACTATGAGCCTTAGGCTTTATTTTCTTATTTACAG
NGCGTCCTCGTTGGAGTGACATCGTCTTTAAACCCCGCGTGGCAATCCCTGACGCACCGCCGTGATGCCCAGGGAAGACAGGGCGACCTGGAAGTCCAACT
NTTCAACAGTGTTTGAACGGAACAGACCCGGGGATTCCCACTGTACTCGCTTCCAGCCGCCTTTACAAGTCTCTCCAGTCGCAGCCTCCGGGACAGATCGG
NTGCAAGCCGGGTATCACAGAGAACCCTCTTTGGACTCTTGGGATTGAAGGTGGATCTTCTCTTCCTGGAAACTCCTTACCTTTCAAATGCAAAGCAAAGA
NCAAAGAGAAGCACTCGCCAAAGACCCCTGGCAAAAAGGCACAACCTCTAGAAGGGCCAGCTGGTCTCAAGGAACACTTTGAAACACCAAACCCCAAAGAT
NTCAAGATGTGAACAGATTCTCATGCCTTCATTTCTCCATCATAGACTCCTCTAGATCCATAAGCTCCCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG
NACTGACCCGGGCAGCCAAGGTGTTGGAGCAGCTCACAGGCCAGACCCCGGTGTTCTCCAAAGCTAGATACACTGTCAGGTCCTTTGGCATCCGGAGAGAT
NGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTTGTTTGGGTTTGGGGCTTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTATCGATCGGTG
NAAAAAGTGATCCCAAAGAAAAATCTAGAGACCGCTCACTGAGCCCTAGGAAAGGAGAAAGTAAAGGTCGGCTCACCATTAAGGCGGGCTCTGGACAAGAT
```

Okay I have no idea how to do this and it feels daunting. Gonna move on and ask Leslie how I could have found it on my own later.

R1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

R2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

Are these actually present? For R1:

```
(QAA) [11:31:15 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0349.talapas.uoregon.edu
$ zcat /projects/bgmp/shared/2017_sequencing/demultiplexed/2_2B_control_S2_L008_R1_001.fastq.gz | grep -A 1 "^@" | grep -v "^@" | grep -v "^--$" | grep -c "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
31917
(QAA) [11:32:17 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0349.talapas.uoregon.edu
$ zcat /projects/bgmp/shared/2017_sequencing/demultiplexed/2_2B_control_S2_L008_R1_001.fastq.gz | grep -A 1 "^@" | grep -v "^@" | grep -v "^--$" | grep -c "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
0
```

This makes sense. Adapter 1 should be present in R1 and Adapter 2 should NOT be present in R1.

For R2:

```
$ zcat /projects/bgmp/shared/2017_sequencing/demultiplexed/2_2B_control_S2_L008_R2_001.fastq.gz | grep -A 1 "^@" | grep -v "^@" | grep -v "^--$" | grep -c "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
0
(QAA) [11:34:04 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0349.talapas.uoregon.edu
$ zcat /projects/bgmp/shared/2017_sequencing/demultiplexed/2_2B_control_S2_L008_R2_001.fastq.gz | grep -A 1 "^@" | grep -v "^@" | grep -v "^--$" | grep -c "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
31965
```

Also makes sense! Okay. Time to use cutadapt.

Documentation:
https://cutadapt.readthedocs.io/en/stable/guide.html

Help:
```
(QAA) [11:34:19 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0349.talapas.uoregon.edu
$ cutadapt -help
cutadapt version 4.9

Copyright (C) 2010 Marcel Martin <marcel.martin@scilifelab.se> and contributors

Cutadapt removes adapter sequences from high-throughput sequencing reads.

Usage:
    cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq

For paired-end reads:
    cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq

Replace "ADAPTER" with the actual sequence of your 3' adapter. IUPAC wildcard
characters are supported. All reads from input.fastq will be written to
output.fastq with the adapter sequence removed. Adapter matching is
error-tolerant. Multiple adapter sequences can be given (use further -a
options), but only the best-matching adapter will be removed.

Input may also be in FASTA format. Compressed input and output is supported and
auto-detected from the file name (.gz, .xz, .bz2). Use the file name '-' for
standard input/output. Without the -o option, output is sent to standard output.

Citation:

Marcel Martin. Cutadapt removes adapter sequences from high-throughput
sequencing reads. EMBnet.Journal, 17(1):10-12, May 2011.
http://dx.doi.org/10.14806/ej.17.1.200

Run "cutadapt --help" to see all command-line options.
See https://cutadapt.readthedocs.io/ for full documentation.

Options:
  -h, --help            Show this help message and exit
  --version             Show version number and exit
  --debug               Print debug log. Use twice to also print DP matrices
  -j CORES, --cores CORES
                        Number of CPU cores to use. Use 0 to auto-detect. Default: 1

Finding adapters:
  Parameters -a, -g, -b specify adapters to be removed from each read (or from R1 if data is paired-
  end. If specified multiple times, only the best matching adapter is trimmed (but see the --times
  option). Use notation 'file:FILE' to read adapter sequences from a FASTA file.

  -a ADAPTER, --adapter ADAPTER
                        Sequence of an adapter ligated to the 3' end (paired data: of the first read).
                        The adapter and subsequent bases are trimmed. If a '$' character is appended
                        ('anchoring'), the adapter is only found if it is a suffix of the read.
  -g ADAPTER, --front ADAPTER
                        Sequence of an adapter ligated to the 5' end (paired data: of the first read).
                        The adapter and any preceding bases are trimmed. Partial matches at the 5' end
                        are allowed. If a '^' character is prepended ('anchoring'), the adapter is only
                        found if it is a prefix of the read.
  -b ADAPTER, --anywhere ADAPTER
                        Sequence of an adapter that may be ligated to the 5' or 3' end (paired data: of
                        the first read). Both types of matches as described under -a and -g are allowed.
                        If the first base of the read is part of the match, the behavior is as with -g,
                        otherwise as with -a. This option is mostly for rescuing failed library
                        preparations - do not use if you know which end your adapter was ligated to!
  -e E, --error-rate E, --errors E
                        Maximum allowed error rate (if 0 <= E < 1), or absolute number of errors for
                        full-length adapter match (if E is an integer >= 1). Error rate = no. of errors
                        divided by length of matching region. Default: 0.1 (10%)
  --no-indels           Allow only mismatches in alignments. Default: allow both mismatches and indels
  -n COUNT, --times COUNT
                        Remove up to COUNT adapters from each read. Default: 1
  -O MINLENGTH, --overlap MINLENGTH
                        Require MINLENGTH overlap between read and adapter for an adapter to be found.
                        Default: 3
  --match-read-wildcards
                        Interpret IUPAC wildcards in reads. Default: False
  -N, --no-match-adapter-wildcards
                        Do not interpret IUPAC wildcards in adapters.
  --action {trim,retain,mask,lowercase,crop,none}
                        What to do if a match was found. trim: trim adapter and up- or downstream
                        sequence; retain: trim, but retain adapter; mask: replace with 'N' characters;
                        lowercase: convert to lowercase; crop: trim up and downstream sequence; none:
                        leave unchanged. Default: trim
  --rc, --revcomp       Check both the read and its reverse complement for adapter matches. If match is
                        on reverse-complemented version, output that one. Default: check only read

Additional read modifications:
  -u LEN, --cut LEN     Remove LEN bases from each read (or R1 if paired; use -U option for R2). If LEN
                        is positive, remove bases from the beginning. If LEN is negative, remove bases
                        from the end. Can be used twice if LENs have different signs. Applied *before*
                        adapter trimming.
  --nextseq-trim 3'CUTOFF
                        NextSeq-specific quality trimming (each read). Trims also dark cycles appearing
                        as high-quality G bases.
  -q [5'CUTOFF,]3'CUTOFF, --quality-cutoff [5'CUTOFF,]3'CUTOFF
                        Trim low-quality bases from 5' and/or 3' ends of each read before adapter
                        removal. Applied to both reads if data is paired. If one value is given, only
                        the 3' end is trimmed. If two comma-separated cutoffs are given, the 5' end is
                        trimmed with the first cutoff, the 3' end with the second.
  --quality-base N      Assume that quality values in FASTQ are encoded as ascii(quality + N). This
                        needs to be set to 64 for some old Illumina FASTQ files. Default: 33
  --poly-a              Trim poly-A tails
  --length LENGTH, -l LENGTH
                        Shorten reads to LENGTH. Positive values remove bases at the end while negative
                        ones remove bases at the beginning. This and the following modifications are
                        applied after adapter trimming.
  --trim-n              Trim N's on ends of reads.
  --length-tag TAG      Search for TAG followed by a decimal number in the description field of the
                        read. Replace the decimal number with the correct length of the trimmed read.
                        For example, use --length-tag 'length=' to correct fields like 'length=123'.
  --strip-suffix STRIP_SUFFIX
                        Remove this suffix from read names if present. Can be given multiple times.
  -x PREFIX, --prefix PREFIX
                        Add this prefix to read names. Use {name} to insert the name of the matching
                        adapter.
  -y SUFFIX, --suffix SUFFIX
                        Add this suffix to read names; can also include {name}
  --rename TEMPLATE     Rename reads using TEMPLATE containing variables such as {id}, {adapter_name}
                        etc. (see documentation)
  --zero-cap, -z        Change negative quality values to zero.

Filtering of processed reads:
  Filters are applied after above read modifications. Paired-end reads are always discarded pairwise
  (see also --pair-filter).

  -m LEN[:LEN2], --minimum-length LEN[:LEN2]
                        Discard reads shorter than LEN. Default: 0
  -M LEN[:LEN2], --maximum-length LEN[:LEN2]
                        Discard reads longer than LEN. Default: no limit
  --max-n COUNT         Discard reads with more than COUNT 'N' bases. If COUNT is a number between 0 and
                        1, it is interpreted as a fraction of the read length.
  --max-expected-errors ERRORS, --max-ee ERRORS
                        Discard reads whose expected number of errors (computed from quality values)
                        exceeds ERRORS.
  --max-average-error-rate ERROR_RATE, --max-aer ERROR_RATE
                        as --max-expected-errors (see above), but divided by length to account for reads
                        of varying length.
  --discard-trimmed, --discard
                        Discard reads that contain an adapter. Use also -O to avoid discarding too many
                        randomly matching reads.
  --discard-untrimmed, --trimmed-only
                        Discard reads that do not contain an adapter.
  --discard-casava      Discard reads that did not pass CASAVA filtering (header has :Y:).

Output:
  --quiet               Print only error messages.
  --report {full,minimal}
                        Which type of report to print: 'full' or 'minimal'. Default: full
  --json FILE           Dump report in JSON format to FILE
  -o FILE, --output FILE
                        Write trimmed reads to FILE. FASTQ or FASTA format is chosen depending on input.
                        Summary report is sent to standard output. Use '{name}' for demultiplexing (see
                        docs). Default: write to standard output
  --fasta               Output FASTA to standard output even on FASTQ input.
  -Z                    Use compression level 1 for gzipped output files (faster, but uses more space)
  --info-file FILE      Write information about each read and its adapter matches into FILE. See the
                        documentation for the file format.
  -r FILE, --rest-file FILE
                        When the adapter matches in the middle of a read, write the rest (after the
                        adapter) to FILE.
  --wildcard-file FILE  When the adapter has N wildcard bases, write adapter bases matching wildcard
                        positions to FILE. (Inaccurate with indels.)
  --too-short-output FILE
                        Write reads that are too short (according to length specified by -m) to FILE.
                        Default: discard reads
  --too-long-output FILE
                        Write reads that are too long (according to length specified by -M) to FILE.
                        Default: discard reads
  --untrimmed-output FILE
                        Write reads that do not contain any adapter to FILE. Default: output to same
                        file as trimmed reads

Paired-end options:
  The -A/-G/-B/-U/-Q options work like their lowercase counterparts, but are applied to R2 (second
  read in pair)

  -A ADAPTER            3' adapter to be removed from R2
  -G ADAPTER            5' adapter to be removed from R2
  -B ADAPTER            5'/3 adapter to be removed from R2
  -U LENGTH             Remove LENGTH bases from R2
  -Q [5'CUTOFF,]3'CUTOFF
                        Quality-trimming cutoff for R2. Default: same as for R1
  -L LENGTH             Shorten R2 to LENGTH. Default: same as for R1
  -p FILE, --paired-output FILE
                        Write R2 to FILE.
  --pair-adapters       Treat adapters given with -a/-A etc. as pairs. Either both or none are removed
                        from each read pair.
  --pair-filter {any,both,first}
                        Which of the reads in a paired-end read have to match the filtering criterion in
                        order for the pair to be filtered. Default: any
  --interleaved         Read and/or write interleaved paired-end reads.
  --untrimmed-paired-output FILE
                        Write second read in a pair to this FILE when no adapter was found. Use with
                        --untrimmed-output. Default: output to same file as trimmed reads
  --too-short-paired-output FILE
                        Write second read in a pair to this file if pair is too short.
  --too-long-paired-output FILE
                        Write second read in a pair to this file if pair is too long.
```

So I want this basic format for paired end reads:

cutadapt -a ADAPT1 -A ADAPT2 -o out1.fastq -p out2.fastq in1.fastq in2.fastq

This outputs a MASSIVE wall of text. The majority of the info I need is near the top but I want to make sure everything is recorded... I think I'm going to output these to two files. cutadaptout_2_2B_control_S2_L008.txt and cutadaptout_19_3F_fox_S14_L008.txt

```
(QAA) [09:28:49 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0349.talapas.uoregon.edu
$ cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o trimmed_2_2B_control_S2_L008_R_001.fastq.gz -p trimmed_2_2B_control_S2_L008_R2_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/2_2B_control_S2_L008_R1_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/2_2B_control_S2_L008_R2_001.fastq.gz > cutadaptout_2_2B_control_S2_L008.txt
Done           00:01:26     5,830,665 reads @  14.9 µs/read;   4.04 M reads/minute
```

```
(QAA) [09:35:06 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0349.talapas.uoregon.edu
$ cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o trimmed_19_3F_fox_S14_L008_R1_001.fastq.gz -p trimmed_19_3F_fox_S14_L008_R2_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/19_3F_fox_S14_L008_R1_001.fastq.gz /projects/bgmp/shared/2017_sequencing/demultiplexed/19_3F_fox_S14_L008_R2_001.fastq.gz > cutadaptout_19_3F_fox_S14_L008.txt
Done           00:03:56    16,348,255 reads @  14.5 µs/read;   4.14 M reads/minute
```

 ------------------
| Wednesday 240825 |
 ------------------

OKAY! I reran the cutadapt for the commands you see above as originally I had a massive wall of text in here. I'm going to use Trimmomatic next.

These are the parameters I need:
    LEADING: quality of 3
    TRAILING: quality of 3
    SLIDING WINDOW: window size of 5 and required quality of 15
    MINLENGTH: 35 bases
    Be sure to output compressed files and clear out any intermediate files.

Some helpful info:

http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

```
(QAA) [09:26:04 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0349.talapas.uoregon.edu
$ trimmomatic -help
Usage: 
       PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
   or: 
       SE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] <inputFile> <outputFile> <trimmer1>...
   or: 
       -version
```

Here's an example command from that helpful link:
java -jar trimmomatic-0.30.jar PE s_1_1_sequence.txt.gz s_1_2_sequence.txt.gz lane1_forward_paired.fq.gz lane1_forward_unpaired.fq.gz lane1_reverse_paired.fq.gz lane1_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

I have paired end, so I need PE, my two sequences, my four output files, then LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35

Running in sbatch, gonna check on it later.

```
TrimmomaticPE: Started with arguments:
 trimmed_2_2B_control_S2_L008_R1_001.fastq.gz trimmed_2_2B_control_S2_L008_R2_001.fastq.gz fw_paired_2_2B_control_S2_L008_R1_001.fq.gz fw_unpaired_2_2B_control_S2_L008_R1_001.fq.gz rv_paired_2_2B_control_S2_L008_R2_001.fq.gz rv_unpaired_2_2B_control_S2_L008_R2_001.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35
Multiple cores found: Using 4 threads
Quality encoding detected as phred33
Input Read Pairs: 5830665 Both Surviving: 5652538 (96.94%) Forward Only Surviving: 133553 (2.29%) Reverse Only Surviving: 4598 (0.08%) Dropped: 39976 (0.69%)
TrimmomaticPE: Completed successfully
	Command being timed: "trimmomatic PE trimmed_2_2B_control_S2_L008_R1_001.fastq.gz trimmed_2_2B_control_S2_L008_R2_001.fastq.gz fw_paired_2_2B_control_S2_L008_R1_001.fq.gz fw_unpaired_2_2B_control_S2_L008_R1_001.fq.gz rv_paired_2_2B_control_S2_L008_R2_001.fq.gz rv_unpaired_2_2B_control_S2_L008_R2_001.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35"
	User time (seconds): 279.27
	System time (seconds): 7.07
	Percent of CPU this job got: 216%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 2:12.45
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 398432
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 44537
	Voluntary context switches: 65418
	Involuntary context switches: 462
	Swaps: 0
	File system inputs: 0
	File system outputs: 864
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

TrimmomaticPE: Started with arguments:
 trimmed_19_3F_fox_S14_L008_R1_001.fastq.gz trimmed_19_3F_fox_S14_L008_R2_001.fastq.gz fw_paired_19_3F_fox_S14_L008_R1_001.fq.gz fw_unpaired_19_3F_fox_S14_L008_R1_001.fq.gz rv_paired_19_3F_fox_S14_L008_R2_001.fq.gz rv_unpaired_19_3F_fox_S14_L008_R2_001.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35
Multiple cores found: Using 4 threads
Quality encoding detected as phred33
Input Read Pairs: 16348255 Both Surviving: 15899268 (97.25%) Forward Only Surviving: 428450 (2.62%) Reverse Only Surviving: 14868 (0.09%) Dropped: 5669 (0.03%)
TrimmomaticPE: Completed successfully
	Command being timed: "trimmomatic PE trimmed_19_3F_fox_S14_L008_R1_001.fastq.gz trimmed_19_3F_fox_S14_L008_R2_001.fastq.gz fw_paired_19_3F_fox_S14_L008_R1_001.fq.gz fw_unpaired_19_3F_fox_S14_L008_R1_001.fq.gz rv_paired_19_3F_fox_S14_L008_R2_001.fq.gz rv_unpaired_19_3F_fox_S14_L008_R2_001.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35"
	User time (seconds): 798.10
	System time (seconds): 19.26
	Percent of CPU this job got: 213%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 6:23.26
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 402092
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 49697
	Voluntary context switches: 186236
	Involuntary context switches: 861
	Swaps: 0
	File system inputs: 0
	File system outputs: 2464
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

```

Now I need to plot these. Installed matplotlib into QAA.

```
(QAA) [03:02:10 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0352.talapas.uoregon.edu
$ conda install matplotlib
```

 ------------------
| Thursday 240829  |
 ------------------

Okay finished plotting stuff. Miiiight redo this in R?

```
(QAA) [03:32:36 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0349.talapas.uoregon.edu
$ ./lengthdistribution_R1R2_plot.py -f1 2_2B_control_S2_L008/fw_paired_2_2B_control_S2_L008_R1_001.fq.gz -f2 2_2B_control_S2_L008/rv_paired_2_2B_control_S2_L008_R2_001.fq.gz -o readlen_2_2B_control_S2_L008_R1R2.png -u True
Getting read lengths...
Making plot...
Done!
```

```
(QAA) [03:34:18 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0349.talapas.uoregon.edu
$ ./lengthdistribution_R1R2_plot.py -f1 19_3F_fox_S14_L008/fw_paired_19_3F_fox_S14_L008_R1_001.fq.gz -f2 19_3F_fox_S14_L008/rv_paired_19_3F_fox_S14_L008_R2_001.fq.gz -o readlen_19_3F_fox_S14_L008_R1R2.png -u True
Getting read lengths...
Making plot...
Done!
```
 ------------------
| Monday 240902  |
 ------------------

## Part 3 - Alignment and strand-specificity

Installing a buncha packages:

```
(QAA) [02:00:20 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0349.talapas.uoregon.edu
$ conda install star numpy matplotlib htseq
```

I need to: "Find publicly available mouse genome fasta files (Ensemble release 112) and generate an alignment database from them. Align the reads to your mouse genomic database using a splice-aware aligner. Use the settings specified in PS8 from Bi621."

Okay. Gonna wget these from Ensembl. Using the "Mus_musculus.GRCm39.dna.primary_assembly.fa.gz" and "Mus_musculus.GRCm39.112.gtf.gz" files.

```
(base) [02:37:09 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0349.talapas.uoregon.edu
$ wget https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
--2024-09-02 14:37:12--  https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
Resolving ftp.ensembl.org (ftp.ensembl.org)... 193.62.193.169
Connecting to ftp.ensembl.org (ftp.ensembl.org)|193.62.193.169|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 806418890 (769M) [application/x-gzip]
Saving to: ‘Mus_musculus.GRCm39.dna.primary_assembly.fa.gz’

Mus_musculus.GRCm39.dna.primary_assem 100%[=========================================================================>] 769.06M   440KB/s    in 30m 33s 

2024-09-02 15:07:46 (430 KB/s) - ‘Mus_musculus.GRCm39.dna.primary_assembly.fa.gz’ saved [806418890/806418890]
```

```
(QAA) [02:47:40 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0352.talapas.uoregon.edu
$ wget https://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz
--2024-09-02 14:53:50--  https://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz
Resolving ftp.ensembl.org (ftp.ensembl.org)... 193.62.193.169
Connecting to ftp.ensembl.org (ftp.ensembl.org)|193.62.193.169|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 32107516 (31M) [application/x-gzip]
Saving to: ‘Mus_musculus.GRCm39.112.gtf.gz’

Mus_musculus.GRCm39.112.gtf.gz        100%[=========================================================================>]  30.62M   448KB/s    in 71s     

2024-09-02 14:55:02 (443 KB/s) - ‘Mus_musculus.GRCm39.112.gtf.gz’ saved [32107516/32107516]
```

Installed samtools as well:

```
(QAA) [02:31:52 /projects/bgmp/adayton/bioinfo] adayton:n0352.talapas.uoregon.edu
$ conda install samtools
```

Made a directory for my star output:

```
(QAA) [02:40:45 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0352.talapas.uoregon.edu
$ mkdir Mus_musculus.GRCm39.dna_ens112_STAR_27.11b
```

Gotta remember star uses all caps. so STAR.

First the database...

```
	Command being timed: "STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ./Mus_musculus.GRCm39.dna_ens112_STAR_27.11b --genomeFastaFiles ./Mus_musculus.GRCm39.dna.primary_assembly.fa --sjdbGTFfile ./Mus_musculus.GRCm39.112.gtf"
	User time (seconds): 5686.37
	System time (seconds): 61.17
	Percent of CPU this job got: 528%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 18:07.66
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 32378200
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 22291086
	Voluntary context switches: 17526
	Involuntary context switches: 5421
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

Now for the alignment...

OKay there's an issue with my fastq file! Some lines are empty!

```
$ zcat trimmed_2_2B_control_S2_L008_R1_001.fastq.gz | grep -A 10 "@K00337:83:HJKJNBBXX:8:1101:1347:2035" 
@K00337:83:HJKJNBBXX:8:1101:1347:2035 1:N:0:CGATCGAT+ATCGATCG

+

@K00337:83:HJKJNBBXX:8:1101:2118:2035 1:N:0:CGATCGAT+ATCGATCG
CATAAACTATCCCTTGTCAGCGTGACACACAATCAGATCTCAGTCCAAATAAAAACAATAGCCAGGTCACAATAAGGCCTAACAGAATAACTA
+
AAF<AFJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJ<7AFJJJJJAJJJJJJJJJJJJJJJFJFJJJJJJJJJJJJFJJJJJJJJJJ<JJJJJ
@K00337:83:HJKJNBBXX:8:1101:2321:2035 1:N:0:CGATCGAT+ATCGATCG
GGGATATGCTCATCACCCCACTTGACCACCAGTGTGTACTCCCCTTTGTCCTTGAGCAGGTAAGAGACACTATAGAGGCGGCTGCCCATGTGTTTCACCAG
+
```

So they trimmed out the whole adapter, which means I need to somehow get these lines out of here. I can get empty lines with "grep "^$""....

How do I make sure I get the whole record? Writing a little python script called clearemptyrecords.py to handle this. I did this command for all 4 files:

```
(base) [05:12:57 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0349.talapas.uoregon.edu
$ ./clearemptyfastqrecords.py -f 19_3F_fox_S14_L008/trimmed_19_3F_fox_S14_L008_R1_001.fastq.gz -o 19_3F_fox_S14_L008/cleaned_19_3F_fox_S14_L00
8_R1_001.fastq.gz -u True
```

Okay that took forever.

Ran staralign.sh on em.

```
	Command being timed: "STAR --runThreadN 8 --runMode alignReads --outFilterMultimapNmax 3 --outSAMunmapped Within KeepPairs --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesCommand zcat --readFilesIn ../2_2B_control_S2_L008/cleaned_2_2B_control_S2_L008_R1_001.fastq.gz ../2_2B_control_S2_L008/cleaned_2_2B_control_S2_L008_R2_001.fastq.gz --genomeDir ../../Mus_musculus.GRCm39.dna_ens112_STAR_27.11b --outFileNamePrefix ../2_2B_control_S2_L008/align_2_2B_control_S2_L008_"
	User time (seconds): 253.04
	System time (seconds): 9.77
	Percent of CPU this job got: 542%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:48.41
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 27321572
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 134699
	Voluntary context switches: 40108
	Involuntary context switches: 762
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
	
	Command being timed: "STAR --runThreadN 8 --runMode alignReads --outFilterMultimapNmax 3 --outSAMunmapped Within KeepPairs --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesCommand zcat --readFilesIn ../19_3F_fox_S14_L008/cleaned_19_3F_fox_S14_L008_R1_001.fastq.gz ../19_3F_fox_S14_L008/cleaned_19_3F_fox_S14_L008_R2_001.fastq.gz --genomeDir ../../Mus_musculus.GRCm39.dna_ens112_STAR_27.11b --outFileNamePrefix ../19_3F_fox_S14_L008/align_2_19_3F_fox_S14_L008_"
	User time (seconds): 834.29
	System time (seconds): 12.40
	Percent of CPU this job got: 686%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 2:03.27
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 27323508
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 211535
	Voluntary context switches: 113244
	Involuntary context switches: 2440
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

Okay I've been working on and off for like 5 hours. Time to enjoy my Monday evening off!

 ------------------
| Tuesday 240903  |
 ------------------

Gonna run my PS8 sammapped script on the aligned files now. Fixed the script to take argparse.

```
(base) [09:13:36 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0349.talapas.uoregon.edu
$ ./sammapped.py -f 2_2B_control_S2_L008/align_2_2B_control_S2_L008_Aligned.out.sam 
Mapped reads:    200241
Unmapped reads:  11388657
(base) [09:13:59 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA] adayton:n0349.talapas.uoregon.edu
$ ./sammapped.py -f 19_3F_fox_S14_L008/align_2_19_3F_fox_S14_L008_Aligned.out.sam 
Mapped reads:    732750
Unmapped reads:  31963176
```

Next is using htseq-count. "Count reads that map to features using htseq-count. You should run htseq-count twice: once with --stranded=yes and again with --stranded=reverse. Use default parameters otherwise."

Wrote a lil slurm script for this as I have to go to class soon.

 ------------------
| Wednesday 240904  |
 ------------------

OK I have to re run this so I can get separate outputs for each file.

```
	Command being timed: "htseq-count --stranded=yes --stranded=reverse /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA/2_2B_control_S2_L008/align_2_2B_control_S2_L008_Aligned.out.sam /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/mouse_data/Mus_musculus.GRCm39.112.gtf"
	User time (seconds): 308.33
	System time (seconds): 2.91
	Percent of CPU this job got: 97%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 5:19.11
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 252556
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 405554
	Voluntary context switches: 1552
	Involuntary context switches: 500
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

```
	Command being timed: "htseq-count --stranded=yes --stranded=reverse /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA/19_3F_fox_S14_L008/align_2_19_3F_fox_S14_L008_Aligned.out.sam /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/mouse_data/Mus_musculus.GRCm39.112.gtf"
	User time (seconds): 887.11
	System time (seconds): 2.77
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 14:50.82
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 253132
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 471103
	Voluntary context switches: 413
	Involuntary context switches: 485
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

```

I did a bunch of git commits for folder cleaning, made an RMarkdown. Time to start analyzing in full!

 ------------------
| Wednesday 240904  |
 ------------------

```
	Command being timed: "htseq-count --stranded=yes /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA/2_2B_control_S2_L008/align_2_2B_control_S2_L008_Aligned.out.sam /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/mouse_data/Mus_musculus.GRCm39.112.gtf"
	User time (seconds): 308.88
	System time (seconds): 2.65
	Percent of CPU this job got: 98%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 5:16.32
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 252572
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 278667
	Voluntary context switches: 1531
	Involuntary context switches: 168
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

		Command being timed: "htseq-count --stranded=reverse /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA/2_2B_control_S2_L008/align_2_2B_control_S2_L008_Aligned.out.sam /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/mouse_data/Mus_musculus.GRCm39.112.gtf"
	User time (seconds): 311.05
	System time (seconds): 2.87
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 5:14.71
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 253296
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 2
	Minor (reclaiming a frame) page faults: 301117
	Voluntary context switches: 418
	Involuntary context switches: 203
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

		Command being timed: "htseq-count --stranded=yes /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA/19_3F_fox_S14_L008/align_2_19_3F_fox_S14_L008_Aligned.out.sam /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/mouse_data/Mus_musculus.GRCm39.112.gtf"
	User time (seconds): 783.21
	System time (seconds): 5.94
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 13:12.76
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 252944
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 395836
	Voluntary context switches: 1000
	Involuntary context switches: 299
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

		Command being timed: "htseq-count --stranded=reverse /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA/19_3F_fox_S14_L008/align_2_19_3F_fox_S14_L008_Aligned.out.sam /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/mouse_data/Mus_musculus.GRCm39.112.gtf"
	User time (seconds): 780.68
	System time (seconds): 6.01
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 13:09.97
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 249288
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 456231
	Voluntary context switches: 443
	Involuntary context switches: 317
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

Went back to ICA4 and used my old command to find the percentage of reads mapped with htseq-count.

First the number of reads that are mapped, unmapped, total # of reads.

```
(base) [01:19:55 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA/htseqcount_data] adayton:n0352.talapas.uoregon.edu
$ awk '$1!~"__" {mapped+=$2} $1~"__" {unmapped+=$2}  END {print mapped, unmapped, (mapped+unmapped)}' 2_2B_control_S2_L008_htseqcount_fwd.txt 
4598 5789851 5794449
(base) [01:19:58 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA/htseqcount_data] adayton:n0352.talapas.uoregon.edu
$ awk '$1!~"__" {mapped+=$2} $1~"__" {unmapped+=$2}  END {print mapped, unmapped, (mapped+unmapped)}' 2_2B_control_S2_L008_htseqcount_rvs.txt 
54800 5739649 5794449
(base) [01:20:28 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA/htseqcount_data] adayton:n0352.talapas.uoregon.edu
$ awk '$1!~"__" {mapped+=$2} $1~"__" {unmapped+=$2}  END {print mapped, unmapped, (mapped+unmapped)}' 19_3F_fox_S14_L008_htseqcount_fwd.txt 
11853 16336110 16347963
(base) [01:20:38 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA/htseqcount_data] adayton:n0352.talapas.uoregon.edu
$ awk '$1!~"__" {mapped+=$2} $1~"__" {unmapped+=$2}  END {print mapped, unmapped, (mapped+unmapped)}' 19_3F_fox_S14_L008_htseqcount_rvs.txt
282927 16065036 16347963
```

Then the percentage of reads:

```
$ awk '$1!~"__" {mapped+=$2} $1~"__" {unmapped+=$2}  END {print (mapped/(mapped+unmapped))*100}' 2_2B_control_S2_L008_htseqcount_fwd.txt 
0.0793518
(base) [01:22:22 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA/htseqcount_data] adayton:n0352.talapas.uoregon.edu
$ awk '$1!~"__" {mapped+=$2} $1~"__" {unmapped+=$2}  END {print (mapped/(mapped+unmapped))*100}' 2_2B_control_S2_L008_htseqcount_rvs.txt 
0.945733
(base) [01:22:25 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA/htseqcount_data] adayton:n0352.talapas.uoregon.edu
$ awk '$1!~"__" {mapped+=$2} $1~"__" {unmapped+=$2}  END {print (mapped/(mapped+unmapped))*100}' 19_3F_fox_S14_L008_htseqcount_fwd.txt 
0.0725044
(base) [01:22:29 /projects/bgmp/adayton/bioinfo/Bi623/Assignments/QAA/QAA/htseqcount_data] adayton:n0352.talapas.uoregon.edu
$ awk '$1!~"__" {mapped+=$2} $1~"__" {unmapped+=$2}  END {print (mapped/(mapped+unmapped))*100}' 19_3F_fox_S14_L008_htseqcount_rvs.txt 
1.73066
```

Looks like this library is strand-specific to the reverse strand.

Okay! I formatted all of this in my R markdown, so I can submit once I figure out the adapter question...
