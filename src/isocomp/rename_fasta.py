
#!/usr/bin/env python3

"""
For some background information:
    Step 0)

    With three input FASTAs (e.g. HG002, HG004, HG005), for each fasta file:

    (0B) Make sure FASTA seq IDs are both unique and have the
    sample name. Proposal: "SampleName_i_existingFastaSeqID"
    whereby i is the increment. <- This script takes care of this part
    of step 0.

    (0C) Align these three FASTAs to the T2T reference, and
    create a BAM file.

Brief description:
    The script takes one FASTA file, rename the sequence
    names to "SampleName_i_existingFastaSeqID" and print the output to
    an output FASTA file (default = <sample>.renamed.fasta, if not specified)

usage:
    step0_rename_fasta_desc.py -i input.fasta -o output.fasta -s sample1
    step0_rename_fasta_desc.py -i input.fasta -s sample1
"""

import argparse
import os

parser = argparse.ArgumentParser(description="Rename FASTA sequence names")

parser.add_argument(
    "-i", "--input", required=True, help="path to input FASTA file", action="store"
)
parser.add_argument(
    "-o", "--output", help="path to output renamed FASTA file", action="store",
)
parser.add_argument(
    "-s",
    "--sample",
    required=True,
    help="sample name to append to sequence name, e.g. sample1",
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
    sample = args.sample
    input_fa = args.input
    if not os.path.exists(input_fa):
        raise SystemError("Error: %s does not exist\n" % input_fa)
    if not args.output:
        output_fa = sample + ".renamed.fasta"
    else:
        output_fa = args.output
    rename_fa_desc(input_fa, output_fa, sample)