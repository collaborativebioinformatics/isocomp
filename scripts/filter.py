import sys
import os
import re

inFile = sys.argv[1] # path to input file
outFile = sys.argv[2] # prefix to output file


for j in range(0,11,1):

    cut = (1000-j)/1000
    out = open(outFile+'%s.tsv' %(str(cut)), 'w')

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

