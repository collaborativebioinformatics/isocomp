"""
    File name: trim_indel_ends.py
    Author: Yutong Qiu
    Date: Nov 20
    Description: For each pair of alinment between reads, if the best alignment includes indels at the beginning and ends
        of the reads, then trim off the indels.
    Output: new set of reads for each sample.
"""

import pandas as pd
from collections import defaultdict, Counter
import re

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

def get_prefix_and_suffix_indels(cigar_str:str) -> tuple:
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
    return prefix_indels, suffix_indels

def trim_seq(prefix_type, prefix_len, suffix_type, suffix_len, seq, num) -> str:
    """Trim input sequence according to suffix and prefix number of indels

    Args:
        prefix_type (str): type of edits at the beginning of the cigar, ("I" or "D")
        prefix_len (int): length of prefix edits
        suffix_type (str):type of edits at the end of the cigar, ("I" or "D")
        suffix_len (int): length of suffix edits
        seq (str): input sequence to be trimmed
        num (int): ranking of the sequence. If a sequence is the first in alignment, deletion means that removing. 
        Otherwise, insertion means trimming.

    Returns:
        str: trimmed sequence
    """
    assert num == 1 or num == 2, "trim_seq(): unknown ranking of input sequence!"

    if num == 1:
        if prefix_type == "D":
            seq = seq[prefix_len:]
        if suffix_type == "D":
            seq = seq[:-suffix_len]
    if num == 1:
        if prefix_type == "I":
            seq = seq[prefix_len:]
        if suffix_type == "I":
            seq = seq[:-suffix_len]
    return seq

def trim_from_cigar(cigar_str, sequence1, sequence2):
    """Given two sequences and the cigar string describing the alignment between them, trim the start and end, if the 
    alignment starts and/or ends with indels.

    Args:
        cigar_str (str): cigar_string
        sequence1(str): aligned sequence 1 to be trimmed
        sequence2(str): aligned sequence 2 to be trimmed
    Returns:
        str: trimmed sequence 1
        str: trimmed sequence 2
    """

    prefix_indels, suffix_indels = get_prefix_and_suffix_indels(cigar_str)

    if len(prefix_indels) == 0 and len(suffix_indels) == 0:
        return sequence1, sequence2
    global total_trim
    total_trim += 1
    # print("Prefix and suffix indels:", prefix_indels, suffix_indels)

    # print("Before trimming:", len(sequence1), len(sequence2))

    prefix_type = ""
    prefix_length = 0
    if len(prefix_indels) != 0:
        for edit in prefix_indels:
            if prefix_type != "":
                assert prefix_type == edit[-1]
            else:
                prefix_type = edit[-1]
            prefix_length += len(edit[:-1])

    suffix_type = ""
    suffix_length = 0
    if len(suffix_indels) != 0:
        for edit in suffix_indels:
            if suffix_type != "":
                assert suffix_type == edit[-1]
            else:
                suffix_type = edit[-1]
            suffix_length += len(edit[:-1])

    sequence1 = trim_seq(prefix_type, prefix_length, suffix_type, suffix_length, sequence1, 1)
    sequence2 = trim_seq(prefix_type, prefix_length, suffix_type, suffix_length, sequence2, 2)

    global correct_trim
    if len(sequence1) == len(sequence2):
        correct_trim += 1

    # print("After trimming:", len(sequence1), len(sequence2))

    return sequence1, sequence2

    

total_trim = 0

input_dir = "/Users/yutongq/SVHackathon2022/"
fname = input_dir+"unique_isoforms_and_aln_stats.tsv"


"""Header
    win_chr
    win_start
    win_end
    total_isoform
    isoform_name
    sample_from
    sample_compared_to
    mapped_start
    isoform_sequence
    selected_alignments
"""
aln_stats_df = pd.read_csv(fname, sep="\t").dropna()

'''
    For each sample:
        for each read in sample:
            will there be multiple reads that are associated with each other?
'''

records = aln_stats_df.to_dict("records")

unique_read_pairs = defaultdict(dict)             # {read_id1: {read_id2: edit distance}}
read_pair_cigar = defaultdict(dict)               # {read_id1: {read_id2: cigar_string}}
read_sequence = dict()                            # {read_id: sequence}

for record in records:
    read_id = record["sample_from"] + "_" + record["isoform_name"]
    compared_to= "_".join(record["selected_alignments"].split("_")[1:3])
    cigar = record["selected_alignments"].split("_")[-1]

    read_sequence[read_id] = record["isoform_sequence"]

    # # check if read_id, compared_to pair is already accounted for.
    # if compared_to in unique_read_pairs and read_id in unique_read_pairs[compared_to]:
    #     continue

    if compared_to not in unique_read_pairs[read_id]:
        read_pair_cigar[read_id][compared_to] = cigar
        unique_read_pairs[read_id].add(compared_to)

multi_pair = Counter()
unique_pair_cigar = defaultdict(dict)

correct_trim = 0

# print(read_pair_cigar)
for read_id in unique_read_pairs:
    multi_pair[len(unique_read_pairs[read_id])] += 1

    if len(unique_read_pairs[read_id]) == 1:
        for compared_to in unique_read_pairs[read_id]:
                cigar = read_pair_cigar[read_id][compared_to]
                seq1 = read_sequence[read_id]
                try:
                    seq2 = read_sequence[compared_to]
                except KeyError:
                    break
                new_seq1, new_seq2 = trim_from_cigar(cigar, seq1, seq2)
                old_ed = ed_from_cigar(read_pair_cigar[read_id][compared_to])
                ## TODO new_ed = 
                break


print("total trim:", total_trim)
print("Correct trim:", correct_trim)


    

    # if len(unique_read_pairs[r]) > 1:
    #     eds = []
    #     for compared_to in unique_read_pairs[r]:
    #         print(r,compared_to+",")
    #         ed = ed_from_cigar(read_pair_cigar[r][compared_to])
    #         eds.append(ed)
    #     print(eds)
print("Multi pair", multi_pair)