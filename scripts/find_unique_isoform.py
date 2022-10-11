# -*- coding: UTF-8 -*-
import os, sys, argparse, re, datetime, pysam
from Bio import SeqIO
#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['inFiles', 'inBed', 'outFile']
optList = ['mapQCut']
sys.stdout.write('\n'); sys.stdout.flush()
# author and version info
usage = sys.argv[0] + ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0 (2022-10-11)'
description = '\nDescription: The program finds unique isoforms by samples and genomic regions'

parser = argparse.ArgumentParser(usage = '\n'.join([usage, author, version, description]), formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, help='string\tinput manifest files.  e.g. /in/samples.csv')
parser.add_argument(posList[1], type=str, help='string\tinput bed file.  e.g. /in/window.bed')
parser.add_argument(posList[2], type=str, help='string\toutput file.  e.g. /out/file.tsv')
# optional arguments
parser.add_argument('-q', '--'+optList[0], type=int, metavar='', default=20, help='integer\tminimum mapQ for reads, default 20')

# strict positional/optional arguments checking
argsDict = vars(parser.parse_args())
sys.stdout.write('\n' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - ' + sys.argv[0] + ' - Parsing Input Arguements...\n\n'); sys.stdout.flush()
for key, value in argsDict.items():
    if key in posList: sys.stdout.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - Required Argument - ' + key +': '+ str(value) + '\n'); sys.stdout.flush()
    if key in optList: sys.stdout.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - Optional Argument - ' + key +': '+ str(value) + '\n'); sys.stdout.flush()
    vars()[key] = value # assign values of arguments into shorthand global variables
sys.stdout.write('\n'); sys.stdout.flush()

#--------------------------------------------------------
# global variables and user-defined functions
#--------------------------------------------------------
# function to config input bed file
def bedConfigFunc(inBed):

    bedDict = {}
    with open(inBed) as f:

        for line in f:
            chr,start,end = line.strip().split()[:3]
            if chr in bedDict: bedDict[chr].append([start,end])
            if not chr in bedDict: bedDict[chr] = [[start,end]]

    return(bedDict)


def readSeqsFunc(inFiles, bedDict):

    with open(inFiles) as f: fields = [ l.strip().split(',') for l in f ]

    outDict = {}
    for chr,beds in bedDict.items():
        for start,end in beds:
            outDict[(chr,start,end)] = {}

    for sample,bam,fa in fields:
        seqDict = {s.id : str(s.seq) for s in SeqIO.parse(fa, 'fasta')}
        samfile = pysam.AlignmentFile(bam, 'rb')

        for chr,beds in bedDict.items():
            for start,end in beds:
                outDict[(chr,start,end)][sample] = {}
                for read in samfile.fetch(chr, int(start), int(end)):
                    qname = read.query_name.split('_')[2]
                    outDict[(chr,start,end)][sample][qname] = {}
                    outDict[(chr,start,end)][sample][qname]['seq'] = seqDict[qname]
                    outDict[(chr,start,end)][sample][qname]['start'] = read.reference_start
                    outDict[(chr,start,end)][sample][qname]['query_length'] = read.query_length
                #outDict[(chr,start,end)][sample]['seqs'] = set([ seqDict[read.query_name.split('_')[2]] for read in samfile.fetch(chr, int(start), int(end)) ])

                #outDict[(chr,start,end)][sample].append(seqDict[name])
            #out.write( '>%s_%s:%s-%s\n%s\n' %(name, chr, start, end, seqDict[name]) )
        samfile.close()

    return outDict


def compareFunc(seqDict, outFile):
    
    out = open(outFile, 'w')

    for locus,allSeqs in seqDict.items():
        for sample,sampleDict in allSeqs.items():
            unique = set([ sampleDict[qname]['seq'] for qname in sampleDict.keys() ])
            #others = []
            #for s in allSeqs.keys():
            #    if s != sample: others += [ [allSeqs[s][q] for q in allSeqs[s].keys()] for s in allSeqs.keys() if s != sample]
            others = [ [allSeqs[s][q]['seq'] for q in allSeqs[s].keys()] for s in allSeqs.keys() if s != sample]

            for other in others: unique = unique.difference(set(other))

            for u in unique:
                for qname,v in sampleDict.items():
                    if v['seq'] == u: out.write('\t'.join(map(str, list(locus) + [qname, u, sample, v['start'],v['query_length']])) +'\n')

    out.close()


#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    bedDict = bedConfigFunc(inBed)
    seqDict = readSeqsFunc(inFiles, bedDict)
    compareFunc(seqDict, outFile)

sys.stdout.write('\n' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - ' + sys.argv[0] + ' - End of Program\n'); sys.stdout.flush()

