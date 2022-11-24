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
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import sys
import argparse
import datetime

#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['manifest','unique_isoform_aln']
optList = ['mode', 'threshold']
sys.stdout.write('\n'); sys.stdout.flush()

# author and version info
usage = sys.argv[0] + ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Yutong Qiu (yutongq@andrew.cmu.edu)'
version = 'Version: 1.0.0.0 (2022-11-23)'
description = '\nDescription: This program trims isoform sequences to reduce false positives. Output trimmed sequences in trimmed.fasta'

parser = argparse.ArgumentParser(usage = '\n'.join([usage, author, version, description]), formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, help='string\tinput manifest csv file in format: sample,fa.  e.g. /in/manifest.csv')
parser.add_argument(posList[1], type=str, help='string\talignment between similar isoforms within each window.  e.g. /in/unique_isoforms_and_aln_stats.tsv')
# optional arguments
parser.add_argument('-m', '--'+optList[0], type=str, metavar='', default='lenient', help='string\tMode of trimming. "Strict" only trims when sequences match exactly after trimming; "Lenient" trims whenever prefix and suffix indel numbers below threshold, defauls to "Lenient"')
parser.add_argument('-t', '--'+optList[1], type=int, metavar='', default=10000, help='int\tThreshold for trimming, defaults to 10000 (trim every prefix and suffix)')

# strict positional/optional arguments checking
argsDict = vars(parser.parse_args())
sys.stdout.write('\n' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - ' + sys.argv[0] + ' - Parsing Input Arguements...\n\n'); sys.stdout.flush()
for key, value in argsDict.items():
    if key in posList: sys.stdout.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - Required Argument - ' + key +': '+ str(value) + '\n'); sys.stdout.flush()
    if key in optList: sys.stdout.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - Optional Argument - ' + key +': '+ str(value) + '\n'); sys.stdout.flush()
    vars()[key] = value # assign values of arguments into shorthand global variables
sys.stdout.write('\n'); sys.stdout.flush()

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

def get_prefix_and_suffix_indels(cigar_str:str, mode="lenient") -> tuple:
    """Get indel edits from the prefix and suffix of the cigar string

    Args:
        cigar_str (str): cigar string

    Returns:
        tuple: (list, list): prefix list of indels and suffix list of indels
    """
    all_edits = extract_edits(cigar_str)
    prefix_indels = []
    suffix_indels = []
    i = 0
    for i in range(len(all_edits)):
        if all_edits[i][-1] != "I" and all_edits[i][-1] != "D":
            break
        prefix_indels.append(all_edits[i])
    j = len(all_edits)-1
    for j in range(len(all_edits)-1, -1, -1):
        if all_edits[j][-1] != "I" and all_edits[j][-1] != "D":
            break
        suffix_indels.append(all_edits[j])
    for idx in range(i,j+1):
        if all_edits[idx][-1] != "=" and mode == "strict":
            return ("", 0), ("",0)

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

def trim_from_cigar(cigar_str, sequence1, sequence2, threshold=10, mode="lenient"):
    """Given two sequences and the cigar string describing the alignment between them, trim the start and end, if the 
    alignment starts and/or ends with indels.

    Args:
        cigar_str (str): cigar_string
        sequence1(str): aligned sequence 1 to be trimmed
        sequence2(str): aligned sequence 2 to be trimmed
        threshold(int): maximum number of bases trimmed
        mode(str): strict or lenient
            - strict: only trim if the trimmed result will be all matches
            - lenient: trim whenever there are indels at the edge of the strings
    Returns:
        str: trimmed sequence 1
        str: trimmed sequence 2
    """

    (prefix_type, prefix_length), (suffix_type, suffix_length) = get_prefix_and_suffix_indels(cigar_str, mode)

    if prefix_length == 0 and suffix_length == 0:
        return sequence1, sequence2

    global total_trim
    total_trim += 1

    sequence1 = trim_seq(prefix_type, prefix_length, suffix_type, suffix_length, sequence1, 1, threshold)
    sequence2 = trim_seq(prefix_type, prefix_length, suffix_type, suffix_length, sequence2, 2, threshold)

    global effective_trim
    if len(sequence1) == len(sequence2) and sequence1 == sequence2:
        effective_trim += 1

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
    count = 0
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        seq_dict[sample_name + "_" + name] = sequence
        count += 1
    print("Read %i sequences in sample %s" % (count, sample_name))
    return seq_dict

# total number of pairs of sequences getting trimmed
total_trim = 0

# total number of pairs of sequences that are identical after getting trimmed
effective_trim = 0

def parseManifest(manifest:str) -> dict:
    """Parse Manifest file

    Args:
        manifest (str): manifest in format sample,fa_path

    Returns:
        dict: read_sequence: {read_id: sequence}
    """
    print("Parsing fasta files...")

    # read manifest
    with open(manifest) as f: fields = [ l.strip().split(',') for l in f ]

    read_sequence = dict()

    for sample,fa in fields:
        read_sequence = read_fasta(fa, sample, read_sequence)
    return read_sequence

def write_fasta(processed_read, read_sequence, mode, threshold):
    records = defaultdict(list)
    for read_id in read_sequence:
        sample_name = read_id.split("_")[0]
        sequence = read_sequence[read_id]
        if read_id in processed_read:
            sequence = processed_read[read_id]
            if sequence == "":
                sequence = read_sequence[read_id]

        record = SeqRecord( Seq(sequence), id=read_id, description="" )
        records[sample_name].append(record)

    for sample_name in records:
        with open("%s_trimmed_%s_%s.fasta" % (sample_name, mode, str(threshold)), "w") as out:
            SeqIO.write(records[sample_name], out, "fasta")
        print("%i sequences in %s" % (len(records[sample_name]), sample_name))
        print("Trimmed written to %s_trimmed_%s_%s.fasta" % (sample_name, mode, str(threshold)))

def main(manifest, unique_isoform_aln, mode, threshold):

    read_sequence = parseManifest(manifest)
    aln_stats_df = pd.read_csv(unique_isoform_aln, sep="\t")
    aln_records = aln_stats_df.to_dict("records")

    unique_read_pairs = defaultdict()
    unique_read_pairs_ed = defaultdict(dict)          # {read_id1: {read_id2: edit distance}}
    read_pair_cigar = defaultdict(dict)               # {read_id1: {read_id2: cigar_string}}

    print("Parsing read comparisons...")
    for record in tqdm(aln_records):
        read_id = record["sample_from"] + "_" + record["isoform_name"]

        unique_read_pairs[read_id] = ""

        try:
            selected_alignments = record["selected_alignments"].split(",")
        except:
            continue
        for alignment in selected_alignments:
            compared_to= alignment.split("__")[1]

            # # # store a pair only once
            # if compared_to in unique_read_pairs_ed and read_id in unique_read_pairs_ed[compared_to]:
            #     continue

            cigar = alignment.split("__")[-1]

            read_pair_cigar[read_id][compared_to] = cigar

            align = edlib.align(read_sequence[compared_to], read_sequence[read_id],mode="NW", task="path")
            ed = align["editDistance"]
            cigar_rev = align["cigar"]

            read_pair_cigar[compared_to][read_id] = cigar_rev

            unique_read_pairs_ed[read_id][compared_to] = ed
            unique_read_pairs_ed[compared_to][read_id] = ed

            # if a read is matched with multiple reads, keep only the most similar pair
            if read_id in unique_read_pairs and unique_read_pairs[read_id] != "":
                opt_compared_to = unique_read_pairs[read_id]
                if ed < unique_read_pairs_ed[read_id][opt_compared_to]:
                    unique_read_pairs[read_id] = compared_to
            else:
                unique_read_pairs[read_id] = compared_to

            if compared_to in unique_read_pairs and unique_read_pairs[compared_to] != "":
                opt_read_id = unique_read_pairs[compared_to]
                if ed < unique_read_pairs_ed[compared_to][opt_read_id]:
                    unique_read_pairs[compared_to] = read_id
            else:
                unique_read_pairs[compared_to] = read_id

    print("Trimming...")

    old_eds = []
    new_eds = []
    length_diff1 = []
    length_diff2 = []

    processed_read = dict()
    for read_id in tqdm(unique_read_pairs):

        compared_to = unique_read_pairs[read_id]

        processed_read[read_id] = ""

        if compared_to == "":
            processed_read[read_id] = read_sequence[read_id]
            continue

        if compared_to in processed_read and read_id in processed_read:
            continue

        processed_read[compared_to] = ""

        cigar = read_pair_cigar[read_id][compared_to]

        seq1 = read_sequence[read_id]
        seq2 = read_sequence[compared_to]
        new_seq1, new_seq2 = trim_from_cigar(cigar, seq1, seq2, threshold=threshold, mode=mode)
        old_ed = ed_from_cigar(read_pair_cigar[read_id][compared_to])
        new_align = edlib.align(new_seq1, new_seq2,mode = "NW", task = "path" )
        new_ed = new_align['editDistance']

        if seq1 != new_seq1 or seq2 != new_seq2:
            old_eds.append(old_ed)
            new_eds.append(new_ed)
            if new_ed == 0:
                length_diff1.append(abs(len(seq1) - len(new_seq1)))
                length_diff1.append(abs(len(seq2) - len(new_seq2)))
            else:
                length_diff2.append(abs(len(seq1) - len(new_seq1)))
                length_diff2.append(abs(len(seq2) - len(new_seq2)))

            # retain the longer sequence
            
            if len(new_seq1) >= len(processed_read[read_id]):
                processed_read[read_id] = new_seq1
            if len(new_seq2) >= len(processed_read[compared_to]):
                processed_read[compared_to] = new_seq2
    
    write_fasta(processed_read, read_sequence, mode, threshold)

    print("Total reads:", len(read_sequence))
    print("Total pairs:", len(unique_read_pairs))
    print("Total trim:", total_trim)
    print("Correct trim:", effective_trim)

    fig, ax = plt.subplots(1)
    ax.scatter(old_eds, new_eds, s=5)
    ax.plot([0, max(old_eds)], [0, max(old_eds)], color="red")
    ax.set_xrange(0, 20)
    ax.set_xlabel("original edit distance")
    ax.set_ylabel("timmed edit distance")
    ax.set_title("Edit distance before and after trimming")
    fig.savefig("ED_trim_%s_%s.png" % (mode, str(threshold)))
    print("Scatterplot saved to ED_trim_%s_%s.png" % (mode, str(threshold)))


    fig,ax = plt.subplots(1)
    ax.hist(length_diff1, range=[0,100],bins=100,label="ed=0")
    ax.hist(length_diff2, range=[0,100],bins=100, alpha=0.5, label="ed!=0")
    ax.set_title("Distribution of trimmed lengths")
    ax.set_xlabel("Number of bases")
    ax.set_ylabel("Number of sequences")
    ax.legend()

    fig.savefig("ED_trim_length_%s_%s.png" % (mode, str(threshold)))
    print("Histogram saved to ED_trim_length_%s_%s.png" % (mode, str(threshold)))

if __name__ == "__main__":
    main(manifest, unique_isoform_aln, mode, threshold)