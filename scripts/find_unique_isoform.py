#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
import datetime
import pysam
import edlib
import numpy
import logging
import logging
from Bio import SeqIO

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %I:%M:%S', \
        level=logging.DEBUG)


logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %I:%M:%S', \
        level=logging.DEBUG)


#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['manifest', 'inBed', 'outPre']
posList = ['manifest', 'inBed', 'outPre']
optList = ['minPercent']
sys.stdout.write('\n'); sys.stdout.flush()
# author and version info
usage = ' <options> ' + ' '.join(posList) + '\n'
usage = ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0 (2022-10-11)'
description = '\nDescription: The program finds unique isoforms by samples and genomic regions'

parser = argparse.ArgumentParser(usage = '\n'.join([usage, author, version, description]), formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, help='string\tinput manifest csv file in format: sample,bam,fa.  e.g. /in/manifest.csv')
parser.add_argument(posList[1], type=str, help='string\tinput bed file.  e.g. /in/window.bed')
parser.add_argument(posList[2], type=str, help='string\toutput file prefix.  e.g. /out/Pre')
parser.add_argument(posList[2], type=str, help='string\toutput file prefix.  e.g. /out/Pre')
# optional arguments
parser.add_argument('-q', '--'+optList[0], type=float, metavar='', default=95, help='float\tminimum percentile for edit distance output (descending order, value must be between 0 and 100 inclusive), default 95')

# strict positional/optional arguments checking
argsDict = vars(parser.parse_args())
logger.info('Parsing Input Arguements...')
logger.info('Parsing Input Arguements...')
for key, value in argsDict.items():
    if key in posList: logger.info('Required Argument - %s: %s' %(key, value))
    if key in optList: logger.info('Optional Argument - %s: %s' %(key, value))
    if key in posList: logger.info('Required Argument - %s: %s' %(key, value))
    if key in optList: logger.info('Optional Argument - %s: %s' %(key, value))
    vars()[key] = value # assign values of arguments into shorthand global variables


#--------------------------------------------------------
# functions
#--------------------------------------------------------

# isoform class
class isoform:

    def __init__(self):
        self.sample = '' # sample of this isoform, e.g. HG002
        self.name = '' # name of this isoform, e.g. PB.6.2
        self.id = '' # unique id of this isoform, e.g. HG002_29_PB.6.2
        self.seq = '' # sequence of this isoform
        self.contig = '' # contig of this isoform aligned to, e.g. NC_060925.1
        self.start = 0 # base 1 coordinate of the alignment start of this isoform


# isoform class
class isoform:

    def __init__(self):
        self.sample = '' # sample of this isoform, e.g. HG002
        self.name = '' # name of this isoform, e.g. PB.6.2
        self.id = '' # unique id of this isoform, e.g. HG002_29_PB.6.2
        self.seq = '' # sequence of this isoform
        self.contig = '' # contig of this isoform aligned to, e.g. NC_060925.1
        self.start = 0 # base 1 coordinate of the alignment start of this isoform


# function to config input bed & bam & fa files
def readSeqsFunc(inBed:str, manifest:str) -> dict[tuple[str,str,str], list[isoform]]:
    """read all input information
def readSeqsFunc(inBed:str, manifest:str) -> dict[tuple[str,str,str], list[isoform]]:
    """read all input information

    Args:
        inBed (str): input bed file in "chr start end ..." format
        manifest (str): input manifest file in "sample_id,path_to_bam,path_to_fasta" format
    Args:
        inBed (str): input bed file in "chr start end ..." format
        manifest (str): input manifest file in "sample_id,path_to_bam,path_to_fasta" format

    Returns:
        dict[tuple[str,str,str], list[isoform]]: configed isoforms by genomic window
    """
    Returns:
        dict[tuple[str,str,str], list[isoform]]: configed isoforms by genomic window
    """

    outDict = {}

    # read windows
    with open(inBed) as f:
        for l in f: outDict[tuple(l.strip().split()[:3])] = []
        for l in f: outDict[tuple(l.strip().split()[:3])] = []

    # read manifest
    with open(manifest) as f: fields = [ l.strip().split(',') for l in f ]

    # read sequences
    for sample,bam,fa in fields:

        # read fa & bam
        seqDict = {s.id : str(s.seq) for s in SeqIO.parse(fa, 'fasta')}
        samfile = pysam.AlignmentFile(bam, 'rb')

        # config sequences
        for locus in outDict.keys():

            checkList = []
            checkList = []
            for read in samfile.fetch(locus[0], int(locus[1]), int(locus[2])):

                #flag = '0'*(12-len(format(read.flag, 'b'))) + format(read.flag, 'b')
                #if read.mapping_quality < 60: continue # filter by mapping quality 60
                #if flag[0] == 1: continue # ignore supplementary alignment
                #if flag[3] == 1: continue # ignore non-primary alignment

                temp = isoform()
                temp.sample = sample
                temp.name = read.query_name.split('_')[2]
                temp.id = read.query_name

                if not temp.id in checkList:
                    checkList.append(temp.id)
                    temp.seq = seqDict[temp.name]
                    temp.contig = locus[0]
                    temp.start = read.reference_start + 1 # convert 0-based to 1-based
                    temp.query_length = read.query_length
                    outDict[locus].append(temp)

        samfile.close()

    return outDict


# function for unique isoform search and N-W ed stats
def compareFunc(seqDict:dict[tuple[str,str,str], list[isoform]], outPre:str, minPercent:float) -> None:
    """find unique isoforms and corresponding N-W alignment stats

    Args:
        seqDict (dict[tuple[str,str,str], list[isoform]]): configed isoform sequences from readSeqsFunc()
        outPre (str): output file prefix for results
        minPercent (float): minimum percentile for edit distance output (descending order)
    """

    outUni = open(outPre+'.unique_isoform.tsv', 'w')
    outDiver = open(outPre+'.isoform_diversity.tsv', 'w')
def compareFunc(seqDict:dict[tuple[str,str,str], list[isoform]], outPre:str, minPercent:float) -> None:
    """find unique isoforms and corresponding N-W alignment stats

    Args:
        seqDict (dict[tuple[str,str,str], list[isoform]]): configed isoform sequences from readSeqsFunc()
        outPre (str): output file prefix for results
        minPercent (float): minimum percentile for edit distance output (descending order)
    """

    outUni = open(outPre+'.unique_isoform.tsv', 'w')
    outDiver = open(outPre+'.isoform_diversity.tsv', 'w')

    '''
    .unique_isoform.tsv
        win_chr                 genomic window chromosome
        win_start               genomic window start
        win_end                 genomic window end
        total_isoform           total number of isoforms (remove duplicates)
        isoform_name            name of the unique isoform
        sample_from             sample this isoform is from
        sample_compared_to      sample this isoform is compared to
        mapped_start            alignment start of this isoform to the reference genome
        isoform_sequence        sequence of this isoform
        selected_alignments     selected pairwise alignments (by percentile edit distance) between this isoform and compared isoforms,
                                in "edit-distance__isoform(compared)-id__isoform(compared)-chromosome__isoform(compared)-mapped-start__isoform(compared)-sequence__pairwise-alignment-cigar"
    .isoform_diversity.tsv
        win_chr                 genomic window chromosome
        win_start               genomic window start
        win_end                 genomic window end
        total_isoform           total number of isoforms (remove duplicates)
        isoform_name            name of the unique isoform
        sample_from             sample this isoform is from
        mapped_start            alignment start of this isoform to the reference genome
        identical_isoforms      other isoforms in this window that are identical to this isoform
        divergent_isoforms      other isoforms in this window that are not identical to this isoform
    .unique_isoform.tsv
        win_chr                 genomic window chromosome
        win_start               genomic window start
        win_end                 genomic window end
        total_isoform           total number of isoforms (remove duplicates)
        isoform_name            name of the unique isoform
        sample_from             sample this isoform is from
        sample_compared_to      sample this isoform is compared to
        mapped_start            alignment start of this isoform to the reference genome
        isoform_sequence        sequence of this isoform
        selected_alignments     selected pairwise alignments (by percentile edit distance) between this isoform and compared isoforms,
                                in "edit-distance__isoform(compared)-id__isoform(compared)-chromosome__isoform(compared)-mapped-start__isoform(compared)-sequence__pairwise-alignment-cigar"
    .isoform_diversity.tsv
        win_chr                 genomic window chromosome
        win_start               genomic window start
        win_end                 genomic window end
        total_isoform           total number of isoforms (remove duplicates)
        isoform_name            name of the unique isoform
        sample_from             sample this isoform is from
        mapped_start            alignment start of this isoform to the reference genome
        identical_isoforms      other isoforms in this window that are identical to this isoform
        divergent_isoforms      other isoforms in this window that are not identical to this isoform
    '''

    outUni.write('\t'.join(['win_chr', \
    outUni.write('\t'.join(['win_chr', \
                         'win_start', \
                         'win_end', \
                         'total_isoform', \
                         'isoform_name', \
                         'sample_from', \
                         'sample_compared_to', \
                         'mapped_start', \
                         'isoform_sequence', \
                         'selected_alignments']) +'\n')

    outDiver.write('\t'.join(['win_chr', \
                         'win_start', \
                         'win_end', \
                         'total_isoform', \
                         'isoform_name', \
                         'sample_from', \
                         'mapped_start', \
                         'identical_isoforms', \
                         'divergent_isoforms']) +'\n')

    for locus,isoforms in seqDict.items():
    outDiver.write('\t'.join(['win_chr', \
                         'win_start', \
                         'win_end', \
                         'total_isoform', \
                         'isoform_name', \
                         'sample_from', \
                         'mapped_start', \
                         'identical_isoforms', \
                         'divergent_isoforms']) +'\n')

    for locus,isoforms in seqDict.items():

        total = len(isoforms)
        samples = sorted(list( set([ s.sample for s in isoforms ]) ))
        total = len(isoforms)
        samples = sorted(list( set([ s.sample for s in isoforms ]) ))

        lookup = {}

        # unique isoforms
        for sample in samples:
            isoformThisSample = [ s.seq for s in isoforms if s.sample == sample ]
        # unique isoforms
        for sample in samples:
            isoformThisSample = [ s.seq for s in isoforms if s.sample == sample ]

            for anotherSample in samples:
                if not anotherSample == sample:
            for anotherSample in samples:
                if not anotherSample == sample:

                    # isoforms in another sample
                    isoformAnotherSample = [ s.seq for s in isoforms if s.sample == anotherSample ]
                    # intersect to obtain unique isoforms in this sample
                    unique = sorted(list( set(isoformThisSample).difference(set(isoformAnotherSample)) ))
                    # isoforms in another sample
                    isoformAnotherSample = [ s.seq for s in isoforms if s.sample == anotherSample ]
                    # intersect to obtain unique isoforms in this sample
                    unique = sorted(list( set(isoformThisSample).difference(set(isoformAnotherSample)) ))

                    # N-W alignment stats
                    for u in unique:
                    # N-W alignment stats
                    for u in unique:

                        # all condidate sequences for N-W alignment
                        condidateSeqs = [ s.seq for s in isoforms if s.seq != u ]
                        condidateContigs = [ s.contig for s in isoforms if s.seq != u ]
                        condidateStarts = [ s.start for s in isoforms if s.seq != u ]
                        condidateNames = [ s.name for s in isoforms if s.seq != u ]
                        # all condidate sequences for N-W alignment
                        condidateSeqs = [ s.seq for s in isoforms if s.seq != u ]
                        condidateContigs = [ s.contig for s in isoforms if s.seq != u ]
                        condidateStarts = [ s.start for s in isoforms if s.seq != u ]
                        condidateNames = [ s.name for s in isoforms if s.seq != u ]

                        # N-W alignment and save for quick lookup
                        for a in condidateSeqs:
                            if (u,a) not in lookup: lookup[(u,a)] = edlib.align(u, a, mode = 'NW', task = 'path')
                        # N-W alignment and save for quick lookup
                        for a in condidateSeqs:
                            if (u,a) not in lookup: lookup[(u,a)] = edlib.align(u, a, mode = 'NW', task = 'path')

                        # alignment stats
                        alns = [ lookup[(u,a)] for a in condidateSeqs ]
                        eds = [ round(aln['editDistance']/len(u),2) for aln in alns ]
                        cigars = [ aln['cigar'] for aln in alns ]
                        # alignment stats
                        alns = [ lookup[(u,a)] for a in condidateSeqs ]
                        eds = [ round(aln['editDistance']/len(u),2) for aln in alns ]
                        cigars = [ aln['cigar'] for aln in alns ]

                        # select alignments by percentile edit distance
                        if alns:
                            cut = numpy.percentile(eds, 100 - minPercent)
                            selectedAlignments = [ '__'.join([str(ed),condidateNames[i],condidateContigs[i],str(condidateStarts[i]),condidateSeqs[i],cigars[i]]) for i,ed in enumerate(eds) if ed <= cut ]
                        else:
                            selectedAlignments = []

                        # output all isoforms in this sample that has sequence "u"
                        for isoform in [ s for s in isoforms if s.sample == sample ]:
                            temp = [total, isoform.name, sample, anotherSample, isoform.start, u, ','.join(selectedAlignments)]
                            if isoform.seq == u: outUni.write('\t'.join( map(str, list(locus)+temp) ) +'\n')


        # all isoform identical/divergent record
        for isoform in isoforms:

            identical,divergent = [],[]

            for s in isoforms:
                if s.id != isoform.id:
                    if s.seq == isoform.seq:
                        identical.append(s.sample+'_'+s.name)
                        # select alignments by percentile edit distance
                        if alns:
                            cut = numpy.percentile(eds, 100 - minPercent)
                            selectedAlignments = [ '__'.join([str(ed),condidateNames[i],condidateContigs[i],str(condidateStarts[i]),condidateSeqs[i],cigars[i]]) for i,ed in enumerate(eds) if ed <= cut ]
                        else:
                            selectedAlignments = []

                        # output all isoforms in this sample that has sequence "u"
                        for isoform in [ s for s in isoforms if s.sample == sample ]:
                            temp = [total, isoform.name, sample, anotherSample, isoform.start, u, ','.join(selectedAlignments)]
                            if isoform.seq == u: outUni.write('\t'.join( map(str, list(locus)+temp) ) +'\n')


        # all isoform identical/divergent record
        for isoform in isoforms:

            identical,divergent = [],[]

            for s in isoforms:
                if s.id != isoform.id:
                    if s.seq == isoform.seq:
                        identical.append(s.sample+'_'+s.name)
                    else:
                        divergent.append(s.sample+'_'+s.name)
                        divergent.append(s.sample+'_'+s.name)

            temp = [total, isoform.name, isoform.sample, isoform.start, ','.join(identical), ','.join(divergent)]
            outDiver.write('\t'.join( map(str, list(locus)+temp) ) +'\n')
            temp = [total, isoform.name, isoform.sample, isoform.start, ','.join(identical), ','.join(divergent)]
            outDiver.write('\t'.join( map(str, list(locus)+temp) ) +'\n')

    outUni.close()
    outDiver.close()
    outUni.close()
    outDiver.close()


#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    seqDict = readSeqsFunc(inBed, manifest)
    compareFunc(seqDict, outPre, minPercent)
    compareFunc(seqDict, outPre, minPercent)

logger.info('End of Program\n')
logger.info('End of Program\n')

