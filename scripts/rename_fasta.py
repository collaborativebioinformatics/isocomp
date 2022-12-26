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
import logging
import sys
import os

logging.getLogger(__name__).addHandler(logging.NullHandler())


__all__ = ['rename_fa_desc','main']

def rename_fa_desc(input_fa: str, output_fa: str, sample: str, i: int = 1) -> None:
    """rename FASTA description to: ">SampleName_i_existingFastaSeqID"
    whereby i is the increment, and print the renamed FASTA entries
    Args:
        input_fa (str): _description_
        output_fa (str): _description_
        sample (str): _description_
        i (int, optional): _description_. Defaults to 1.
    """
    # TODO docstring

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


def parse_args(args=None):
    Description = "Rename FASTA sequence names" # TODO write more thorough description
    Epilog = "USAGE: isocomp rename_fasta -i genome.fasta -o renamed_genome.fasta -s sample1"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
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
    return parser.parse_args(args)

def main(args=None):
    
    logging.debug('cmd ling arguments: {args}')
    args = parse_args(args)

    # Check inputs
    logging.info('checking input...')
    input_path_list = [args.input]
    for input_path in input_path_list:
        if not os.path.exists(input_path):
            raise FileNotFoundError("Input file DNE: %s" %input_path)

    sample = args.sample
    input_fa = args.input
    
    if not args.output:
        output_fa = sample + ".renamed.fasta"
    else:
        output_fa = args.output
    
    rename_fa_desc(input_fa, output_fa, sample)

if __name__ == "__main__":
    sys.exit(main())