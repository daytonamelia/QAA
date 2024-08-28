#!/usr/bin/env python 

import bioinfo
import argparse
import gzip
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(description="Makes histogram of 'Mean Quality Score vs Base Pair' for given file and read length.")
    parser.add_argument("-f", "--file", help="Name of input file", type=str, required=True)
    parser.add_argument("-l", "--readlen", help="Length of reads in file", type=int, required=True)
    parser.add_argument("-o", "--outputfile", help="Name of output histogram file", type=str, default="perbasehist.png")
    parser.add_argument("-y", "--ylim", help="Maximum y limit for plot - used for plot fiddling.", type=int, default=40)
    parser.add_argument("-u", "--unzip", help="Pass as True if input file needs to be unzipped.", type=bool, default=False)
    return parser.parse_args()

def mean_phred(file: str, readlen: int, compress: bool=False) -> list:
    """Returns a list of mean phred cores per base for a file given the length of reads in the file."""
    # initialize list
    phred_list = []
    for i in range(0,readlen):
        phred_list.append(0.0)
    # sum phred scores at each base
    if compress:
        with gzip.open(file, 'rt') as fh:
            line_counter = 0
            print("{file} open!")
            for line in fh:
                line = line.strip()
                if line_counter%4 == 3:
                    for i, nt in enumerate(line):
                        phred_list[i] += bioinfo.convert_phred(nt)
                if line_counter%1000000 == 0:
                    print(f"Line: {line_counter}")
                line_counter += 1
    else:
        with open(file, 'r') as fh:
            line_counter = 0
            for line in fh:
                line = line.strip()
                if line_counter%4 == 3:
                    for i, nt in enumerate(line):
                        phred_list[i] += bioinfo.convert_phred(nt)
                if line_counter%1000000 == 0:
                    print(f"Line: {line_counter}")
                line_counter += 1
    # returns mean sums
    return [phred_sum/(line_counter/4) for phred_sum in phred_list]

if __name__ == "__main__":
    args = get_args()
    print("Finding mean phreds...")
    phred_list = mean_phred(args.file, args.readlen, args.unzip)
    print("Making plot...")
    fig, ax = plt.subplots()
    ax.plot(range(0,args.readlen), phred_list)
    ax.set(ylim=(0,args.ylim))
    ax.set_xticks(range(0,args.readlen+1), minor = True)
    ax.set_title('Mean Quality Score vs Base Pair')
    ax.set_xlabel('Base Pair (bp)')
    ax.set_ylabel('Mean Quality Score')
    plt.savefig(args.outputfile)
    print("Done!")
