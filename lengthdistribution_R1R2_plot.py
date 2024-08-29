#!/usr/bin/env python 

import argparse
import gzip
import matplotlib.pyplot as plt
import numpy as np
import time

def get_args():
    parser = argparse.ArgumentParser(description="Makes plot of read length distributions for R1 and R2.")
    parser.add_argument("-f1", "--fileR1", help="Name of input file for R1", type=str, required=True)
    parser.add_argument("-f2", "--fileR2", help="Name of input file for R2", type=str, required=True)
    parser.add_argument("-o", "--outputfile", help="Name of output histogram file", type=str, default="lendist_R1R2.png")
    parser.add_argument("-u", "--unzip", help="Pass as True if input file needs to be unzipped.", type=bool, default=False)
    return parser.parse_args()

def readlen_count(file:str, unzip: bool = False) -> dict:
    """Returns a dictionary of the read length counts for a fastq file. Use unzip if it is a gzipped file."""
    linecount = 0
    read_lengths = {}
    if unzip:
        with gzip.open(file, 'rt') as fh:
            for line in fh:
                line = line.strip()
                if linecount % 4 == 1:
                    readlen = len(line)
                    if readlen not in read_lengths:
                        read_lengths[readlen] = 1
                    else:
                        read_lengths[readlen] += 1
                linecount += 1
    else:
        with open(file, 'rt') as fh:
            for line in fh:
                line = line.strip()
                if linecount % 4 == 1:
                    readlen = len(line)
                    if readlen not in read_lengths:
                        read_lengths[readlen] = 1
                    else:
                        read_lengths[readlen] += 1
                linecount += 1
    return read_lengths

args = get_args()
print("Getting read lengths...", flush = True)
R1_counts = readlen_count(args.fileR1, args.unzip)
R2_counts = readlen_count(args.fileR2, args.unzip)
longest_R1 = max(R1_counts.keys())
longest_R2 = max(R2_counts.keys())
longest_read = max(longest_R1, longest_R2)

print("Making plot...", flush = True)
fig, ax = plt.subplots()
bar_width = 0.25
ax.bar(R1_counts.keys(), R1_counts.values(), width=bar_width, label = "R1 Reads")
ax.bar(R2_counts.keys(), R2_counts.values(), width=bar_width, label = "R2 Reads")
ax.set_xticks(range(0,longest_read+1), minor = True)
ax.set_title('Read Length Distribution')
ax.set_xlabel('Read Length')
plt.savefig(args.outputfile)
print("Done!")
