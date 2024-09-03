#!/usr/bin/env python 

import argparse
import gzip


def get_args():
    parser = argparse.ArgumentParser(description="Clears empty records in a fastq file.")
    parser.add_argument("-f", "--file", help="Name of input file", type=str, required=True)
    parser.add_argument("-o", "--output", help="Name of output histogram file", type=str, required=True)
    parser.add_argument("-u", "--unzip", help="Pass as True if input file needs to be unzipped.", type=bool, default=False)
    return parser.parse_args()

args = get_args()
if args.unzip:
    with gzip.open(args.file, 'rt') as rf, gzip.open(args.output, "wt") as wf:
        record = []
        for line in rf:
            line = line.strip()
            if line.startswith("@"):
                if "" not in record:
                    for recordline in record:
                        wf.write(f"{recordline}\n")
                record = []
                record.append(line)
            else:
                record.append(line)
        if "" not in recordline:
            for recordline in record:
                wf.write(f"{recordline}\n")
else:
    with open(args.file, 'r') as rf, open(args.output, "w") as wf:
        record = []
        for line in rf:
            line = line.strip()
            if line.startswith("@"):
                for recordline in record:
                    wf.write(f"{recordline}\n")
                record = []
                record.append(line)
            else:
                record.append(line)
        for recordline in record:
                    wf.write(f"{recordline}\n")


