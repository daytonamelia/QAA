#!/usr/bin/env python

import re
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="Counts mapped and unmapped reads in a file and prints them to stdout.")
    parser.add_argument("-f", "--file", help="Name of input file", type=str, required=True)
    return parser.parse_args()

args = get_args()
mapped = 0
unmapped = 0
with open(args.file, "r") as fh:
    for line in fh:
        line = line.strip("\n")
        if line.startswith("@"):
            continue
        linefields = re.split("\t", line)
        flag = int(linefields[1])
        if((flag & 4) != 4) and ((flag & 256) != 256):
            mapped += 1
        else:
            if ((flag & 256) != 256):
                unmapped += 1

print("Mapped reads:\t", mapped)
print("Unmapped reads:\t", unmapped)