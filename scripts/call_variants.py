#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import os
import sys
import argparse
import datetime
import re

#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['inFile', 'outFile']
optList = ['N']
sys.stdout.write('\n'); sys.stdout.flush()
# author and version info
usage = sys.argv[0] + ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0 (2022-10-30)'
description = '\nDescription: The program calls SNPs and InDels from the unique isoform file'

parser = argparse.ArgumentParser(usage = '\n'.join([usage, author, version, description]), formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, help='string\tinput file.  e.g. /in/unique_isoforms.tsv')
parser.add_argument(posList[1], type=str, help='string\toutput file.  e.g. /out/unique_isoforms_variants.tsv')
# optional arguments
parser.add_argument('-n', '--'+optList[0], type=int, metavar='', default=10, help='float\tmaximum bases for indel combination (inclusive), default 10')

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

# cigar class
class cigar:

    def __init__(self):
        self.string = '' # cigar string
        self.ref = '' # ref sequence
        self.query = '' # query sequence
        self.refID = '' # ref id
        self.queryID = '' # query_id
        self.parsed = [] # parsed cigar
        self.SNPs = [] # called SNPs
        self.InDels = [] # called InDels


    # function to call SNPs and InDels
    def call(self):

        self.parsed = re.findall(r'(\d+)(\w)', self.string.replace('=','M'))
        self.parsed = [ [int(c[0]),c[1]] for c in self.parsed ]

        refPos, queryPos = 1, 1

        for c in self.parsed:

            # matches
            if c[1] == 'M':
                refPos += c[0]
                queryPos += c[0]
            # mismatches
            elif c[1] == 'X':
                self.SNPs.append([ 'snp', refPos, refPos+c[0], queryPos, queryPos+c[0], self.ref[refPos-1:refPos+c[0]-1], self.query[queryPos-1:queryPos+c[0]-1] ])
                refPos += c[0]
                queryPos += c[0]
            # deletions in query
            elif c[1] == 'D':
                self.InDels.append([ 'del', refPos, refPos+c[0], queryPos, queryPos, self.ref[refPos-1:refPos+c[0]-1], '.' ])
                refPos += c[0]
            # insertions in query
            elif c[1] == 'I':
                self.InDels.append([ 'ins', refPos, refPos, queryPos, queryPos+c[0], '.', self.query[queryPos-1:queryPos+c[0]-1] ])
                queryPos += c[0]
            else:
                sys.exit('unknown cigar notation')


    # function to merge close InDels and remove condensed SNPs
    def refine(self, cut):

        updateSNPs, updateInDels, last = [], [], []

        # combine close ins or del
        for current in self.InDels:

            # handle the starting case
            if last == []:
                updateInDels.append(current)
                last = updateInDels[-1]
                continue

            # retrieve the last InDel
            last = updateInDels[-1]
            typeL, refStartL, refEndL, queryStartL, queryEndL, refL, altL = last
            typeC, refStartC, refEndC, queryStartC, queryEndC, refC, altC = current

            # current InDel to be merged with the last one
            if (typeL == typeC or typeL == 'delins') and refStartC-refEndL <= cut and queryStartC-queryEndL <= cut:
                updateInDels[-1] = [ 'delins', refStartL, refEndC, queryStartL, queryEndC, self.ref[refStartL-1:refEndC-1], self.ref[queryStartL-1:queryEndC-1] ]
            # current InDel not to be merged
            else:
                updateInDels.append(current)

        # check SNPs
        for snp in self.SNPs:
            for indel in updateInDels:

                # SNPs contained in a merged delins
                typeS, refStartS, refEndS, queryStartS, queryEndS, refS, altS = snp
                typeI, refStartI, refEndI, queryStartI, queryEndI, refI, altI = indel

                if refStartS >= refStartI and refEndS <= refEndI and queryStartS >= queryStartI and queryEndS <= queryEndI:
                    updateSNPs.append(snp)
                    break

        self.SNPs = [ s for s in self.SNPs if s not in updateSNPs ]
        self.InDels = updateInDels


def main():

    out = open(outFile, 'w')
    #out.write('\t'.join(['#ref_id','query_id','type','ref','alt','ref_start','ref_end','query_start','query_end','ref_seq','query_seq']) +'\n')
    out.write('\t'.join(['#ref_id','query_id','type','ref','alt','ref_pos']) +'\n')

    with open(inFile) as f:

        used = [('','')]*100

        for index,line in enumerate(f):

            if index == 0: continue

            # parse input, skip entries with total==1
            try:
                chr, start, end, total, query, sampleFrom, sampleTo, mapStart, querySeq, alns = line.strip().split('\t')
                alns = alns.split(',')
            except:
                continue

            # skip duplicated entries
            if (sampleFrom, query) in used: continue
            used = used[1:]
            used.append((sampleFrom, query))

            # obtain all candidate isoform alignments
            cigars = [ a.split('__')[-1] for a in alns ]
            refSeqs = [ a.split('__')[-2] for a in alns ]
            refIDs = [ a.split('__')[-3] for a in alns ]

            # parse cigar and call variants
            for i,c in enumerate(cigars):

                target = cigar()
                target.string = c
                target.ref = refSeqs[i]
                target.query = querySeq
                target.refID = refIDs[i]
                target.queryID = sampleFrom+'_'+query

                # call variants
                target.call()
                # refine indels
                target.refine(N)

                for v in target.SNPs + target.InDels:
                    #out.write('\t'.join(map(str, [target.refID,target.queryID,v[0],v[-2],v[-1],v[1],v[2],v[3],v[4],target.ref,target.query])) +'\n')
                    out.write('\t'.join(map(str, [target.refID,target.queryID,v[0],v[-2],v[-1],v[1]])) +'\n')

#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    main()


sys.stdout.write('\n' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - ' + sys.argv[0] + ' - End of Program\n'); sys.stdout.flush()

