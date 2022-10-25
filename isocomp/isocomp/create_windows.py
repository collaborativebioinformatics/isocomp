#!/usr/bin/env python3
'''
Author: Yutong Qiu
Date: Oct 11, 2022

Requires: 
- 100bp windows with coverage for each sample
    expected fields: ["chr", "start", "end", "per base coverage"]
- transcript/RNA boundaries
    expected fields: ["chr", "start", "end", "geneName", "0", "-", "geneID", "biotype"]
'''
import pybedtools as pybed
import pandas as pd
from collections import defaultdict
import sys, os
import errno
import warnings
warnings.filterwarnings("ignore")

def print_help():
    """Print help
    """
    print("USAGE:\tpython %s <bp_cov_names> <gene_window_name>" % sys.argv[0])
    print()
    print("\tbp_cov_names \t\tcomma-delimited paths to 100bp windows with coverage for each sample")
    print("\t\t\t\t   expected fields: [chr, start, end, per base coverage]")
    print("\tgene_window_name \tpath to transcript/RNA boundaries")
    print("\t\t\t\t   expected fields: [chr, start, end, geneName, 0, -, geneID, biotype]")


def get_all_windows(gene_df:pd.DataFrame, bp_df:pd.DataFrame) -> pd.DataFrame:
    """From gene boundaries and 100 bp nonzero coverage windows, produce a merged window df

    Args:
        gene_df (pd.DataFrame): one window per gene, > 0.05 avg coverage
        bp_df (pd.DataFrame): one window per 100 bp, > 0.05 avg coverage

    Returns:
        pd.DataFrame: merged windows df
    """
    gene_df_pointer = 0
    coverage_pointer = 0
    
    all_windows = defaultdict(list)

    bp_chr, bp_start, bp_end = 0,0,0
    
    while gene_df_pointer < len(gene_df.index):
        curr_chr, curr_start, curr_end = 0,0,0
        
        if coverage_pointer < len(bp_df.index):
            bp_chr = bp_df.loc[coverage_pointer]["chr"]
            bp_start = bp_df.loc[coverage_pointer]["start"]
            bp_end = bp_df.loc[coverage_pointer]["end"]

        gene_chr = gene_df.loc[gene_df_pointer]["chr"]
        gene_start = gene_df.loc[gene_df_pointer]["start"]
        gene_end = gene_df.loc[gene_df_pointer]["end"]
        
        if bp_chr == gene_chr:
            curr_chr = bp_chr
            curr_start = min(gene_start, bp_start)
            
            # if exhausted 100bp windows, use gene
            if coverage_pointer >= len(bp_df.index)-1:
                curr_start = gene_start

            if curr_start == gene_start or gene_start - bp_end < 100:
                if bp_end - gene_end > 100 or gene_start - bp_end < 100:
                    curr_end = gene_end
                else:
                    curr_end = bp_end
                    coverage_pointer += 1
                    
                all_windows["chr"].append(curr_chr)
                all_windows["start"].append(curr_start)
                all_windows["end"].append(curr_end)
                gene_df_pointer += 1


            else:
                curr_end = bp_end
                while(bp_start <= gene_start and bp_end - curr_end < 200 and bp_chr == curr_chr and coverage_pointer < len(bp_df.index)-1):
                    
                    curr_end = bp_end
                    
                    coverage_pointer += 1
                    
                    bp_chr = bp_df.loc[coverage_pointer]["chr"]
                    bp_start = bp_df.loc[coverage_pointer]["start"]
                    bp_end = bp_df.loc[coverage_pointer]["end"]

                
                all_windows["chr"].append(curr_chr)
                all_windows["start"].append(curr_start)
                all_windows["end"].append(curr_end)
                
        # in the case where chromosomes are different:
        else:
            if gene_chr < bp_chr:
                all_windows["chr"].append(gene_chr)
                all_windows["start"].append(gene_start)
                all_windows["end"].append(gene_end)
                gene_df_pointer += 1
            else:                
                all_windows["chr"].append(bp_chr)
                all_windows["start"].append(bp_start)
                all_windows["end"].append(bp_end)
                coverage_pointer += 1
                
                
    # deal with remaining 100bp windows
    while (coverage_pointer < len(bp_df.index)-1):
        curr_chr = bp_chr
        curr_start = bp_start
        curr_end = bp_end
        while(bp_end - curr_end < 200 and bp_chr == curr_chr and coverage_pointer < len(bp_df.index)-1):

            curr_end = bp_end

            coverage_pointer += 1

            bp_chr = bp_df.loc[coverage_pointer]["chr"]
            bp_start = bp_df.loc[coverage_pointer]["start"]
            bp_end = bp_df.loc[coverage_pointer]["end"]
            
        if coverage_pointer == len(bp_df.index) - 1:
            curr_end = bp_end

        all_windows["chr"].append(curr_chr)
        all_windows["start"].append(curr_start)
        all_windows["end"].append(curr_end)
            
    all_windows = pd.DataFrame(all_windows)
    return all_windows


def further_merge(input_windows: pd.DataFrame) -> pd.DataFrame:
    """Further merge windows if they are < 100 bp apart

    Args:
        input_windows (pd.DataFrame): a set of windows from merging gene boundaries and 100bp windows with non zero coverage

    Returns:
        pd.DataFrame: Merged windows
    """
    curr_idx = 0

    new_windows = defaultdict(list)

    curr_chr = input_windows.loc[curr_idx]["chr"]
    curr_start = input_windows.loc[curr_idx]["start"]
    curr_end = input_windows.loc[curr_idx]["end"]
    curr_idx += 1
    while curr_idx < len(input_windows.index):

        add = True

        next_chr = input_windows.loc[curr_idx]["chr"]
        next_start = input_windows.loc[curr_idx]["start"]
        next_end = input_windows.loc[curr_idx]["end"]

        if next_chr == curr_chr:
            if next_start - curr_end < 100:
                curr_end = max(next_end, curr_end)
                add = False

        if add:
            new_windows["chr"].append(curr_chr)
            new_windows["start"].append(curr_start)
            new_windows["end"].append(curr_end)

            curr_chr = next_chr
            curr_start = next_start
            curr_end = next_end
        curr_idx += 1
    new_windows = pd.DataFrame(new_windows)
    return new_windows

if __name__ == "__main__":
        
    if len(sys.argv) != 3:
        print_help()
        exit(1)

    bp_cov_names = sys.argv[1].split(",")
    gene_window_name = sys.argv[2]

    assert len(bp_cov_names) == 3

    int_dir = "intermediate_files/"
    out_dir = "output_files/"

    try:
        os.mkdir(int_dir)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass

    try:
        os.mkdir(out_dir)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass

    # merge 100bp coverage into one file
    print("Merging 100bp coverage files...")
    command = "paste %s | cut -f 1,2,3,4,8,12 > " % " ".join(bp_cov_names) +int_dir+"HG_merged_100bp_coverage.bed" 
    print(command)
    os.system(command)

    # filter out 100bp windows with < 0.05 sum average coverage
    print("Filtering out < 0.05 coverage bins...")

    coverage_100bp_df = pd.read_csv(int_dir+"HG_merged_100bp_coverage.bed", 
                                    sep="\t", header=None, 
                                    names=["chr", "start", "end", "hg002", "hg004", "hg005"])
    coverage_100bp_df["sum_cov"] = coverage_100bp_df["hg002"] + coverage_100bp_df["hg004"] + coverage_100bp_df["hg005"]
    coverage_100bp_df = coverage_100bp_df[coverage_100bp_df["sum_cov"] > 0.05]
    coverage_100bp_df.to_csv(int_dir+"merged_100bp_coverage.nonzero.bed", sep="\t", header=None, index=False)


    # obtain coverage within each gene region
    print("Obtaining coverage per gene...")
    coverage_100bp_bed = pybed.BedTool(int_dir+"merged_100bp_coverage.nonzero.bed")
    gene_window_bed = pybed.BedTool(gene_window_name)
    gene_coverage_df = coverage_100bp_bed.intersect(gene_window_bed, wa=True, wb=True).to_dataframe()
    gene_coverage_df.columns = ["chr1", "start1", "end1", "hg002", "hg004", "hg005","sum_cov", "chr", "start", "end", "geneName", "0", "-", "geneID", "biotype"]

    new_gene_coverage_df = pd.DataFrame()
    new_gene_coverage_df["gene"] = gene_coverage_df["geneName"].drop_duplicates()
    new_gene_coverage_df = new_gene_coverage_df.set_index("gene")
    new_gene_coverage_df["sum_cov"] = gene_coverage_df.groupby("geneName").sum()["sum_cov"]
    new_gene_coverage_df["chr"] = gene_coverage_df.groupby("geneName").min()["chr"]
    new_gene_coverage_df["start"] = gene_coverage_df.groupby("geneName").min()["start"]
    new_gene_coverage_df["end"] = gene_coverage_df.groupby("geneName").max()["end"]
    new_gene_coverage_df = new_gene_coverage_df.reset_index()
    new_gene_coverage_df = new_gene_coverage_df[["chr", "start", "end", "sum_cov", "gene"]]

    new_gene_coverage_df.to_csv(out_dir+"gene_overlap_sum_coverage.bed",sep="\t", header=False, index=False)

    # get nonzero 100bp regions that does not overlap with gene region
    print("Getting 100bp windows outside of gene regions...")
    gene_coverage_bed= pybed.BedTool(out_dir+"gene_overlap_sum_coverage.bed")
    no_overlap_100bp_df = coverage_100bp_bed.intersect(gene_coverage_bed, v=True, wa=True).to_dataframe()
    no_overlap_100bp_df.columns = ["chr","start", "end","hg002","hg004","hg005","sum_cov"]

    no_overlap_100bp_df.to_csv(out_dir+"no_overlap_100bp.nonzero.bed", sep="\t", header=False, index=False)

    # merge gene windows and 100bp windows that do not overlap with genes
    print("Creating windows...")
    merged_windows = get_all_windows(gene_coverage_df, no_overlap_100bp_df)
    # new_all.to_csv("all_output.bed", sep="\t", header=None, index=None)
    print("Final merging...")
    merged_windows = further_merge(merged_windows)
    merged_windows.to_csv(out_dir+"4_merged_windows.bed", sep="\t", header=None, index=None)

    print("Output in "+out_dir+"4_merged_windows.bed")