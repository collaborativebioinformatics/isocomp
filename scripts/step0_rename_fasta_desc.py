#!/usr/bin/env python3

import argparse
import os

parser = argparse.ArgumentParser(description="Rename FASTA sequence names")

parser.add_argument(
    "-i", "--input", required=True, help="path to input FASTA file", action="store"
)
parser.add_argument(
    "-o",
    "--output",
    required=True,
    help="path to output renamed FASTA file",
    action="store",
)
parser.add_argument(
    "-s",
    "--sample",
    required=True,
    help="sample name to append to sequence name",
    action="store",
)


def rename_fa_desc(input_fa: str, output_fa: str, sample: str, i: int = 1):
    # rename FASTA description to: ">SampleName_i_existingFastaSeqID"
    # whereby i is the increment, and print the renamed FASTA entries
    with open(input_fa, "r") as f, open(output_fa, "w") as r:
        print("Reading: ", sample)
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                desc = ">" + sample + "_" + str(i) + "_" + line[1:]
                i += 1
                r.write(desc + "\n")
            else:
                r.write(line + "\n")
        print("Total no of entries for", sample, " : ", i - 1)
        print("-" * 30)


if __name__ == "__main__":
    args = parser.parse_args()
    input_fa = args.input
    if not os.path.exists(input_fa):
        raise SystemError("Error: %s does not exist\n" % input_fa)
    output_fa = args.output
    sample = args.sample
    rename_fa_desc(input_fa, output_fa, sample)
