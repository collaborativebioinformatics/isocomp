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
posList = ['inFile', 'outPre']
optList = []
sys.stdout.write('\n'); sys.stdout.flush()
# author and version info
usage = sys.argv[0] + ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0 (2022-10-31)'
description = '\nDescription: The program filters unique isoform outputs'

parser = argparse.ArgumentParser(usage = '\n'.join([usage, author, version, description]), formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, help='string\tinput file.  e.g. /in/unique_isoforms.tsv')
parser.add_argument(posList[1], type=str, help='string\toutput prefix.  e.g. /out/unique_isoforms_filter')
# optional arguments

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


#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":


    for j in range(0,11,1):

        cut = (1000-j)/1000
        out = open(outPre+'%s.tsv' %(str(cut)), 'w')

        with open(inFile) as f:

            for i,line in enumerate(f):

                # header line
                if i == 0:
                    out.write(line)
                    continue

                # length of this unique isoform
                length = len(line.strip().split()[8])

                # get all selected alignments
                selected = line.strip().split()[-1].split(',')

                # get all cigars (also replacing '=' by 'M' for easier re matching)
                cigars = [ s.split('_')[-1].replace('=','M') for s in selected ]

                # parse all cigars
                cigars = [ re.findall(r'(\d+)(\w)', c) for c in cigars ]

                # get number of matches in the cigars
                matches = [ sum([ int(a[0]) for a in c if a[1]=='M' ]) for c in cigars ]

                # get the best match
                percent = max([ round(m/length,5) for m in matches ])

                # filter
                if percent <= cut: out.write(line)

        out.close()


sys.stdout.write('\n' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - ' + sys.argv[0] + ' - End of Program\n'); sys.stdout.flush()

