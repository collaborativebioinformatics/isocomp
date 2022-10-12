# -*- coding: UTF-8 -*-
import os
import sys
import argparse
import datetime
import pysam
import edlib
import numpy
from Bio import SeqIO

#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['manifest', 'inBed', 'outFile']
optList = ['minPercent']
sys.stdout.write('\n'); sys.stdout.flush()
# author and version info
usage = sys.argv[0] + ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0 (2022-10-11)'
description = '\nDescription: The program finds unique isoforms by samples and genomic regions'

parser = argparse.ArgumentParser(usage = '\n'.join([usage, author, version, description]), formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, help='string\tinput manifest csv file in format: sample,bam,fa.  e.g. /in/manifest.csv')
parser.add_argument(posList[1], type=str, help='string\tinput bed file.  e.g. /in/window.bed')
parser.add_argument(posList[2], type=str, help='string\toutput file.  e.g. /out/file.tsv')
# optional arguments
parser.add_argument('-q', '--'+optList[0], type=float, metavar='', default=95, help='float\tminimum percentile for edit distance output (descending order, value must be between 0 and 100 inclusive), default 95')

# strict positional/optional arguments checking
argsDict = vars(parser.parse_args())
sys.stdout.write('\n' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - ' + sys.argv[0] + ' - Parsing Input Arguements...\n\n'); sys.stdout.flush()
for key, value in argsDict.items():
    if key in posList: sys.stdout.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - Required Argument - ' + key +': '+ str(value) + '\n'); sys.stdout.flush()
    if key in optList: sys.stdout.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - Optional Argument - ' + key +': '+ str(value) + '\n'); sys.stdout.flush()
    vars()[key] = value # assign values of arguments into shorthand global variables
sys.stdout.write('\n'); sys.stdout.flush()


#--------------------------------------------------------
# functions
#--------------------------------------------------------

# function to config input bed & bam & fa files
def readSeqsFunc(inBed, manifest):

    '''
    read all input information.

    @parameter inBed - input bed file in "chr start end ..." format.
    @parameter manifest - input manifest file in "sample_id,path_to_bam,path_to_fasta" format.

    '''

    outDict = {}
    '''
    {
    locus1: {
        sample1: {
            qname1: {
                seq:
                start:
                length:
            }
            qname2: {...}
        }
        sample2: {...}
    locus2: {...}
    ...
    }
    '''

    # read windows
    with open(inBed) as f:
        for l in f: outDict[tuple(l.strip().split()[:3])] = {}

    # read manifest
    with open(manifest) as f: fields = [ l.strip().split(',') for l in f ]

    # read sequences
    for sample,bam,fa in fields:

        # read fa & bam
        seqDict = {s.id : str(s.seq) for s in SeqIO.parse(fa, 'fasta')}
        samfile = pysam.AlignmentFile(bam, 'rb')

        # config sequences
        for locus in outDict.keys():

            outDict[locus][sample] = {}

            for read in samfile.fetch(locus[0], int(locus[1]), int(locus[2])):
                qname = read.query_name.split('_')[2]
                outDict[locus][sample][qname] = {}
                outDict[locus][sample][qname]['seq'] = seqDict[qname]
                outDict[locus][sample][qname]['start'] = read.reference_start
                outDict[locus][sample][qname]['query_length'] = read.query_length

        samfile.close()

    return outDict


# function for unique isoform search and N-W ed stats
def compareFunc(seqDict, outFile, minPercent):

    '''
    find unique isoforms and corresponding N-W alignment stats.

    @parameter seqDict - configed sequences from readSeqsFunc().
    @parameter outFile - output file for results.
    @parameter minPercent - minimum percentile for edit distance output (descending order).

    '''

    out = open(outFile, 'w')
    out.write('\t'.join(['win_chr', \
                         'win_start', \
                         'win_end', \
                         'total_isoform', \
                         'isoform_name', \
                         'sample_from', \
                         'sample_compared_to', \
                         'mapped_start', \
                         'isoform_sequence', \
                         'selected_alignments']) +'\n')

    for locus,allSamples in seqDict.items():

        # all isoforms of the locus
        '''allSeqs, allNames = [], []
        for s,d in allSamples.items():
            allSeqs += [ d[q]['seq'] for q in d.keys() ]
            allNames += [ s+'_'+q for q in d.keys() ]
        total = len(allSeqs)
        '''
        total = sum([ len([ allSamples[sample][qname]['seq'] for qname in allSamples[sample].keys() ]) for sample in allSamples.keys() ])

        lookup = {}

        for sample,sampleDict in allSamples.items():

            # isoforms in this sample
            this = [ sampleDict[q]['seq'] for q in sampleDict.keys() ]

            #lookup = {}

            for anotherSample in allSamples.keys():

                if anotherSample == sample: continue

                # isoforms in another sample
                another = [ allSamples[anotherSample][q]['seq'] for q in allSamples[anotherSample].keys() ]
                # intersect to obtain unique isoforms in this sample
                unique = set(this).difference(set(another))

                # N-W alignment stats
                for u in unique:

                    # all condidate sequences for N-W alignment
                    allSeqs, allNames = [], []
                    for s,d in allSamples.items():
                        allSeqs += [ d[q]['seq'] for q in d.keys() if d[q]['seq'] != u]
                        allNames += [ s+'_'+q for q in d.keys() if d[q]['seq'] != u ]

                    # N-W alignment and save for quick lookup
                    for a in allSeqs:
                        if (u,a) not in lookup: lookup[(u,a)] = edlib.align(u, a, mode = 'NW', task = 'path')

                    # alignment stats
                    alns = [ lookup[(u,a)] for a in allSeqs ]
                    eds = [ round(aln['editDistance']/len(u),2) for aln in alns ]
                    cigars = [ aln['cigar'] for aln in alns ]

                    # select alignments by percentile edit distance
                    if alns:
                        cut = numpy.percentile(eds, 100 - minPercent)
                        picks = [ str(ed)+'_'+allNames[i]+'_'+cigars[i] for i,ed in enumerate(eds) if ed <= cut ]
                    else:
                        picks = []

                    # output
                    for qname,v in sampleDict.items():
                        temp = [total, qname, sample, anotherSample, v['start'], u, ','.join(picks)]
                        if v['seq'] == u: out.write('\t'.join( map(str, list(locus)+temp) ) +'\n')

    out.close()


#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    seqDict = readSeqsFunc(inBed, manifest)
    compareFunc(seqDict, outFile, minPercent)

sys.stdout.write('\n' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - ' + sys.argv[0] + ' - End of Program\n'); sys.stdout.flush()

