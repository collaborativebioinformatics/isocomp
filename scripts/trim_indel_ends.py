"""
    File name: trim_indel_ends.py
    Author: Yutong Qiu
    Date: Nov 20
    Description: For each pair of alinment between reads, if the best alignment includes indels at the beginning and ends
        of the reads, then trim off the indels.
    Input: 
        - old sets of reads for each sample
        - aln_stats
    Output: 
        - new set of reads for each sample.

    Trimming logic alternatives:
        - Trim prefix and suffix regardless of the rest of the alignment --- set a threshold
        - Trim prefix and suffix only when the rest of the alignment is all matches --- set a threshold
"""

import pandas as pd
from collections import defaultdict, Counter
import re
import edlib
from Bio import SeqIO
import matplotlib.pyplot as plt
from tqdm import tqdm

def extract_edits(cigar_str) -> list:
    """Given a cigar string, return separated edits with counts

    Args:
        cigar_str (str):

    Returns:
        list: a list of tuples of (count, edit_type)
    """
    all_edits = [(int(s[:-1]), s[-1]) for s in re.findall(r'\d+[=,I,D,X]', cigar_str)]
    return all_edits

def ed_from_cigar(cigar_str) -> int:
    """Compute edit distance from cigar strings

    Args:
        cigar_str (str): cigar string

    Returns:
        int: edit distance that is equal to total length of the string minus the nu
    """
    all_edits = extract_edits(cigar_str)
    total_length = sum([e[0] for e in all_edits])
    match_count = sum([e[0] for e in all_edits if e[1] == "="])
    return total_length - match_count

def get_len_type(indel_list:list) -> tuple[str, str]:
    """Given a list of indels, return the total length and type

    Args:
        indel_list (list): list of indels

    Returns:
        tuple[str, str]: type, length
    """

    indel_type = ""
    indel_length = 0
    if len(indel_list) != 0:
        for edit in indel_list:
            if indel_type != "":
                assert indel_type == edit[-1]
            else:
                indel_type = edit[-1]
            indel_length += edit[0]
    return indel_type, indel_length

def get_prefix_and_suffix_indels(cigar_str:str) -> tuple:
    """Get indel edits from the prefix and suffix of the cigar string

    Args:
        cigar_str (str): cigar string

    Returns:
        tuple: (list, list): prefix list of indels and suffix list of indels
    """
    all_edits = extract_edits(cigar_str)
    prefix_indels = []
    suffix_indels = []
    for i in range(len(all_edits)):
        if all_edits[i][-1] != "I" and all_edits[i][-1] != "D":
            break
        prefix_indels.append(all_edits[i])
    for i in range(len(all_edits)-1, -1, -1):
        if all_edits[i][-1] != "I" and all_edits[i][-1] != "D":
            break
        suffix_indels.append(all_edits[i])

    # return prefix_indels, suffix_indels[::-1]
    return get_len_type(prefix_indels), get_len_type(suffix_indels[::-1])

def trim_seq(prefix_type, prefix_len, suffix_type, suffix_len, seq, num, threshold=10) -> str:
    """Trim input sequence according to suffix and prefix number of indels

    Args:
        prefix_type (str): type of edits at the beginning of the cigar, ("I" or "D")
        prefix_len (int): length of prefix edits
        suffix_type (str):type of edits at the end of the cigar, ("I" or "D")
        suffix_len (int): length of suffix edits
        seq (str): input sequence to be trimmed
        num (int): ranking of the sequence. If a sequence is the first in alignment, deletion means that removing. 
        Otherwise, insertion means trimming.
        threshold (int): maximum number of bases trimmed

    Returns:
        str: trimmed sequence
    """
    assert num == 1 or num == 2, "trim_seq(): unknown ranking of input sequence!"

    if num == 2:
        if prefix_type == "D":
            seq = seq[min(threshold, prefix_len):]
        if suffix_type == "D":
            seq = seq[:-min(threshold, suffix_len)]
    if num == 1:
        if prefix_type == "I":
            seq = seq[min(threshold, prefix_len):]
        if suffix_type == "I":
            seq = seq[:-min(threshold, suffix_len)]
    return seq

def trim_from_cigar(cigar_str, sequence1, sequence2, threshold=10):
    """Given two sequences and the cigar string describing the alignment between them, trim the start and end, if the 
    alignment starts and/or ends with indels.

    Args:
        cigar_str (str): cigar_string
        sequence1(str): aligned sequence 1 to be trimmed
        sequence2(str): aligned sequence 2 to be trimmed
        threshold(int): maximum number of bases trimmed
    Returns:
        str: trimmed sequence 1
        str: trimmed sequence 2
    """

    (prefix_type, prefix_length), (suffix_type, suffix_length) = get_prefix_and_suffix_indels(cigar_str)

    if prefix_length == 0 and suffix_length == 0:
        return sequence1, sequence2

    global total_trim
    total_trim += 1

    sequence1 = trim_seq(prefix_type, prefix_length, suffix_type, suffix_length, sequence1, 1, threshold)
    sequence2 = trim_seq(prefix_type, prefix_length, suffix_type, suffix_length, sequence2, 2, threshold)

    global correct_trim
    if len(sequence1) == len(sequence2) and sequence1 == sequence2:
        correct_trim += 1

    return sequence1, sequence2

    
def read_fasta(fname:str, sample_name:str, seq_dict:dict) -> dict:
    """Parse fasta file. add sequences to input sequence dict

    Args:
        fname (str): bam file name
        sample_name(str): sample name
        seq_dict(dict): input and output dict

    Returns:
        dict: dictionary matching read name to sequence
    """
    fasta_sequences = SeqIO.parse(open(fname),'fasta')
    for fasta in tqdm(fasta_sequences):
        name, sequence = fasta.id, str(fasta.seq)
        seq_dict[sample_name + "_" + name] = sequence
    return seq_dict

total_trim = 0
correct_trim = 0

def main():
    input_dir_aln = "/Users/yutongq/SVHackathon2022/isocomp/data/"
    input_dir_fasta = "/Users/yutongq/SVHackathon2022/"

    aln_fname = input_dir_aln+"unique_isoforms_and_aln_stats.tsv"
    fasta_1_fname = input_dir_fasta+"HG002_MM2_corrected.fasta"
    fasta_2_fname = input_dir_fasta+"HG004_MM2_corrected.fasta"
    fasta_3_fname = input_dir_fasta+"HG005_MM2_corrected.fasta"

    read_sequence = dict()

    print("Parsing fasta files...")
    read_sequence = read_fasta(fasta_1_fname, "HG002", read_sequence)
    read_sequence = read_fasta(fasta_2_fname, "HG004", read_sequence)
    read_sequence = read_fasta(fasta_3_fname, "HG005", read_sequence)

    aln_stats_df = pd.read_csv(aln_fname, sep="\t").dropna()
    aln_records = aln_stats_df.to_dict("records")

    unique_read_pairs = defaultdict()
    unique_read_pairs_ed = defaultdict(dict)          # {read_id1: {read_id2: edit distance}}
    read_pair_cigar = defaultdict(dict)               # {read_id1: {read_id2: cigar_string}}

    print("Parsing read comparisons...")
    for record in tqdm(aln_records):
        read_id = record["sample_from"] + "_" + record["isoform_name"]
        selected_alignments = record["selected_alignments"].split(",")
        for alignment in selected_alignments:
            compared_to= alignment.split("__")[1]

            # # store a pair only once
            if compared_to in unique_read_pairs_ed and read_id in unique_read_pairs_ed[compared_to]:
                continue

            cigar = alignment.split("_")[-1]

            read_pair_cigar[read_id][compared_to] = cigar
            ed = edlib.align(read_sequence[read_id], read_sequence[compared_to])["editDistance"]
            unique_read_pairs_ed[read_id][compared_to] = ed

            # if a read is matched with multiple reads, keep only the most similar pair
            if read_id in unique_read_pairs:
                opt_compared_to = unique_read_pairs[read_id]
                if ed < unique_read_pairs_ed[read_id][opt_compared_to]:
                    unique_read_pairs[read_id] = compared_to
            else:
                unique_read_pairs[read_id] = compared_to

    print("Trimming...")

    old_eds = []
    new_eds = []
    length_diff1 = []
    length_diff2 = []

    processed_read = set()
    for read_id in tqdm(unique_read_pairs):
        compared_to = unique_read_pairs[read_id]
        cigar = read_pair_cigar[read_id][compared_to]

        seq1 = read_sequence[read_id]
        seq2 = read_sequence[compared_to]
        new_seq1, new_seq2 = trim_from_cigar(cigar, seq1, seq2, threshold=100000)
        old_ed = ed_from_cigar(read_pair_cigar[read_id][compared_to])
        new_align = edlib.align(new_seq1, new_seq2,mode = "NW", task = "path" )
        new_ed = new_align['editDistance']
        new_cigar = new_align["cigar"]

        if seq1 != new_seq1 or seq2 != new_seq2:
            old_eds.append(old_ed)
            new_eds.append(new_ed)
            if new_ed == 0:
                length_diff1.append(abs(len(seq1) - len(new_seq1)))
                length_diff1.append(abs(len(seq2) - len(new_seq2)))
            else:
                length_diff2.append(abs(len(seq1) - len(new_seq1)))
                length_diff2.append(abs(len(seq2) - len(new_seq2)))
        processed_read.add(read_id)
        processed_read.add(compared_to)

    
    # for read_id in read_sequence:
    #     if read_id not in processed_read:
    #         print(read_id)

        # print("Length:", len(seq1), len(seq2), len(new_seq1), len(new_seq2))

    print("Total reads:", len(read_sequence))
    print("Total pairs:", len(unique_read_pairs))
    print("Total trim:", total_trim)
    print("Correct trim:", correct_trim)

    fig, ax = plt.subplots(1)
    ax.scatter(old_eds, new_eds, s=5)
    ax.plot([0, max(old_eds)], [0, max(old_eds)], color="red")
    ax.set_xlabel("original edit distance")
    ax.set_ylabel("timmed edit distance")
    ax.set_title("Edit distance before and after trimming")
    fig.savefig("ED_trim.png")


    fig,ax = plt.subplots(1)
    ax.hist(length_diff1, range=[0,100],bins=100)
    ax.hist(length_diff2, range=[0,100],bins=100, alpha=0.5)

    fig.savefig("ED_trim_length.png")

if __name__ == "__main__":
    main()

    # seq1="GTCTTTCCTCATTTACAGCTAGCTCTGCCAGTCTGTGTCTTTGCTGTGTCTGTGGATGTGAGTGTCCCTGCTGGTCTGTAGGAGATGGTATTTTGGGGGCAGCTGCAAGGGAAAAAACTGATCTGCTACTTACACAGCGCCGTAGCCTCAGCCTGAGGGTCTTCAGGTTCTCCCCCAGGGAGTTCACATGCGCCTTGATGTCTGGGTCTTGGTTCTCAGCTTGGGGCATCACCTCCTCCAGGTAAAACTGGATCATCTCAGACAAGGCTTGGCAACCCAGGTAACCCTAAGGGCAGGAGCCAAAGAGCGGTGCTTGCACACACTGACAGGAGTCCAAGAATGTGCACTGAGGGAGCGTTTCCGCACAGATCTGCGTGTTCCTTACCACTCACACATGTGCACACACATATCCATGTGTGTGTGCCAGTGCTTTGGGGCTCTGTTCCACGGGGTGCAAACCCTAGACGAAAAAAGACTGTAATTCAACTCTTCAGTAGCCAGGATTGCCTGTCTCCAACCCCTTCTGAATCTAGTTATGTTGGATACTGGGAGAGAAGCATGTAGACTCTTTGGGCCTGACCTTTCTGATTTTTCTCTGGAAGGGCTTCACTGAAGGGGCACAGAAGAATAATGTTTCTACAGCATAGATCCAGGGGCTCCTCCCCACCAGGAACCCAGCTTAGCAGGCAGGCTGTTCATGTGGAACTCGCCAAAACCCCGTGCAGGTGGATGTGCGGAAACACCAAAGAACTGCACGATTGTCCCTGCCTTTGGAAGGTCTCATTCCCTGTCATTCCTTCATTCAGTCAACAAACTTGGAGGACCCAGTATGTGCTACATTTTGAGGATAAAGTGATGAACAACACAGACAAGATCCCTGCCCTTGTGGAGTCTATATTCAAACCAGAGCTTGACCTTCAAACTGAGAGACCAACTTAAAAATGGGGCCATCTGGCCACTAACTGTTATACACCCTCTGCTTGGGACACACGGTCTCAGCCCTCCTCCACCTTCAAGTCTGATAGGAGAGAGATTTCAGAGTGGCCCAAACAATCTAAAGTAAAAGTGAACAGAGTTCAAAATAAGTCCTCAAGGGCTGAGGGGAGACCAAGGGGAAGATAATGCGGTGAGGGCACTGAGAGAGGGATTCCTAGAGGAATTGTACCTTGCTTTGAATTTTAGCCCTAAGCTGCAAATCCTAGCTAACTTCTAAGTAGATGCTGACAGCTGGCTTCTATTTGGTCTGCCATATTGAATGGCTTTGTCCCCAGGACATGGTGGAGATTTAGACAACCCAGGGCAGAAAATCCTCTAGGAAAGGAAGTCTGAACTACACCGCTGGCTTTCTGAAAGGACAGGATTTGAGGAGCTGGGATGGTAGTAGGGGAGAGGGCAGCTGGGTGGGGTTTCTTCCAGGAAACAATGTGGATGCAAGATGCAGTCCTCGTCCCTGGCCTGCCCTTACCGCTTTAGTATTTCAGATGAGTCAGCTCTCACTCTGCCACCTGCTTTCTGTAACTCTGCAGCAGGCCAGCAGCCTCCAGTCCACCCTCTCCTCCTCCCAGGGCTCTGGGGTGAAGCCCTCCCCAGGAACTCCCTGGCTCCCAGTACCCATGGGAGAAGCTCTTTTCAACTGACATGATTTATTTTAACAATAAGAGGAAAAAACCACATCCGATTTCCACAGAATTAGAGTTCCTCCTGACACGCCTGGACTTGCCAGTTGGGAAGGACACTCTCAAGCCATCCTCTATGCCATTTGCCGTTCTGCTGTTCCTTTGCTTGGAGGGAAGAGATGCCTGTTTCTATAGGGAGGAAGTTTGAAGTGAGAAGAGAGGACAAAGATTCTGTGAATGGAAAGAGCACTGGACTGAGAGTTCAGGGATGTGGCCTCCGGTCCCAGCTCCACCACTACTCAGTTGTCTGACGAAGGAAACAGAAAGTTGTGATACAATGTGGGTGACACATTGCCTAGAAAGAAGTGCCGTTGTGATACGCTTATGTTGGTGGGATGGGGGAAAGAAATAACAGGTTTGTATGGAGTGTTATGAAAGAATTAAATCTTAACCTTCCATTGGGGTAAGCTGATGGAAACAACCATAGCAAAGGACAGACTGAATATTTTTTTATTCTCTTTATAGAAAATAACAGTAAAAAAAAGTATTGACATAGGAGAAGGAAACAAAAAGTTGTCAGGAGTTAATGAATAGAAATATTATTTTTTTCTGGATTTTATGGTGTTTGTGTTATTTGTTAGCTTTTATGAATTTGTAATTTGTTGTAATTTCTTTCTAAATAAGTATTTACTTTTGTCTCT"

    # seq2 = "TAGCTCTGCCAGTCTGTGTCTTTGCTGTGTCTGTGGATGTGAGTGTCCCTGCTGGTCTGTAGGAGATGGTATTTTGGGGGCAGCTGCAAGGGAAAAAACTGATCTGCTACTTACACAGCGCCGTAGCCTCAGCCTGAGGGTCTTCAGGTTCTCCCCCAGGGAGTTCACATGCGCCTTGATGTCTGGGTCTTGGTTCTCAGCTTGGGGCATCACCTCCTCCAGGTAAAACTGGATCATCTCAGACAAGGCTTGGCAACCCAGGTAACCCTAAGGGCAGGAGCCAAAGAGCGGTGCTTGCACACACTGACAGGAGTCCAAGAATGTGCACTGAGGGAGCGTTTCCGCACAGATCTGCGTGTTCCTTACCACTCACACATGTGCACACACATATCCATGTGTGTGTGCCAGTGCTTTGGGGCTCTGTTCCACGGGGTGCAAACCCTAGACGAAAAAAGACTGTAATTCAACTCTTCAGTAGCCAGGATTGCCTGTCTCCAACCCCTTCTGAATCTAGTTATGTTGGATACTGGGAGAGAAGCATGTAGACTCTTTGGGCCTGACCTTTCTGATTTTTCTCTGGAAGGGCTTCACTGAAGGGGCACAGAAGAATAATGTTTCTACAGCATAGATCCAGGGGCTCCTCCCCACCAGGAACCCAGCTTAGCAGGCAGGCTGTTCATGTGGAACTCGCCAAAACCCCGTGCAGGTGGATGTGCGGAAACACCAAAGAACTGCACGATTGTCCCTGCCTTTGGAAGGTCTCATTCCCTGTCATTCCTTCATTCAGTCAACAAACTTGGAGGACCCAGTATGTGCTACATTTTGAGGATAAAGTGATGAACAACACAGACAAGATCCCTGCCCTTGTGGAGTCTATATTCAAACCAGAGCTTGACCTTCAAACTGAGAGACCAACTTAAAAATGGGGCCATCTGGCCACTAACTGTTATACACCCTCTGCTTGGGACACACGGTCTCAGCCCTCCTCCACCTTCAAGTCTGATAGGAGAGAGATTTCAGAGTGGCCCAAACAATCTAAAGTAAAAGTGAACAGAGTTCAAAATAAGTCCTCAAGGGCTGAGGGGAGACCAAGGGGAAGATAATGCGGTGAGGGCACTGAGAGAGGGATTCCTAGAGGAATTGTACCTTGCTTTGAATTTTAGCCCTAAGCTGCAAATCCTAGCTAACTTCTAAGTAGATGCTGACAGCTGGCTTCTATTTGGTCTGCCATATTGAATGGCTTTGTCCCCAGGACATGGTGGAGATTTAGACAACCCAGGGCAGAAAATCCTCTAGGAAAGGAAGTCTGAACTACACCGCTGGCTTTCTGAAAGGACAGGATTTGAGGAGCTGGGATGGTAGTAGGGGAGAGGGCAGCTGGGTGGGGTTTCTTCCAGGAAACAATGTGGATGCAAGATGCAGTCCTCGTCCCTGGCCTGCCCTTACCGCTTTAGTATTTCAGATGAGTCAGCTCTCACTCTGCCACCTGCTTTCTGTAACTCTGCAGCAGGCCAGCAGCCTCCAGTCCACCCTCTCCTCCTCCCAGGGCTCTGGGGTGAAGCCCTCCCCAGGAACTCCCTGGCTCCCAGTACCCATGGGAGAAGCTCTTTTCAACTGACATGATTTATTTTAACAATAAGAGGAAAAAACCACATCCGATTTCCACAGAATTAGAGTTCCTCCTGACACGCCTGGACTTGCCAGTTGGGAAGGACACTCTCAAGCCATCCTCTATGCCATTTGCCGTTCTGCTGTTCCTTTGCTTGGAGGGAAGAGATGCCTGTTTCTATAGGGAGGAAGTTTGAAGTGAGAAGAGAGGACAAAGATTCTGTGAATGGAAAGAGCACTGGACTGAGAGTTCAGGGATGTGGCCTCCGGTCCCAGCTCCACCACTACTCAGTTGTCTGACGAAGGAAACAGAAAGTTGTGATACAATGTGGGTGACACATTGCCTAGAAAGAAGTGCCGTTGTGATACGCTTATGTTGGTGGGATGGGGGAAAGAAATAACAGGTTTGTATGGAGTGTTATGAAAGAATTAAATCTTAACCTTCCATTGGGGTAAGCTGATGGAAACAACCATAGCAAAGGACAGACTGAATATTTTTTTATTCTCTTTATAGAAAATAACAGTAAAAAAAAGTATTGACATAGGAGAAGGAAACAAAAAGTTGTCAGGAGTTAATGAATAGAAATATTATTTTTTTCTGGATTTTATGGTGTTTGTGTTATTTGTTAGCTTTTATGAATTTGTAATTTGTTGTAATTTCTTTCTAAATAAGTATTTACTTTTGTCTCTAATTTTGTCCATCTATTTTTGAATTCAATTTTCAAATTCAAAACGCTTCTCTCGGCCGGGCACGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGTCGAGGAGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACAGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCAGGCGAGGTGGTGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCATGAACCCCGGGGGGCGGAGCCTGCAGTGAGCCGAGATCGTGCCACTGCACTCCAGCCTGGGTGACAGCGAGACTCCGTCTC"

    # cigar = "1I1=8I1=6I3=2I2=2I2291=333D"

    # new_seq1, new_seq2 = trim_from_cigar(cigar, seq1, seq2, 1000)

    # print(len(seq1), len(seq2))
    # print(len(new_seq1), len(new_seq2))