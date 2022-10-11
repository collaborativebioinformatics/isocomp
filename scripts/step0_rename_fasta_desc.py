#!/usr/bin/env python

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


samples = ["HG002", "HG004", "HG005"]
for sample in samples:
    input_fa = sample + "_corrected.fasta"
    output_fa = sample + "_corrected.renamed.fasta"
    rename_fa_desc(input_fa, output_fa, sample)
