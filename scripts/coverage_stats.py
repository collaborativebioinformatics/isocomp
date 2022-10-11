# -*- coding: UTF-8 -*-
import os, sys, argparse, re, datetime
#--------------------------------------------------------
# parser
#--------------------------------------------------------
# positional/optional arguments
posList = ['inBam', 'inBed', 'samtools', 'outBed']
optList = ['mapQCut']
sys.stdout.write('\n'); sys.stdout.flush()
# author and version info
usage = sys.argv[0] + ' <options> ' + ' '.join(posList) + '\n'
author = 'Author: Bida Gu (bida.beta.gu@gmail.com)'
version = 'Version: 1.0.0.0 (2022-10-11)'
description = '\nDescription: The program generate bam coverage statistics for given bed regions'

parser = argparse.ArgumentParser(usage = '\n'.join([usage, author, version, description]), formatter_class=argparse.RawTextHelpFormatter)
# positional arguments
parser.add_argument(posList[0], type=str, help='string\tinput bam file.  e.g. /in/sample.bam')
parser.add_argument(posList[1], type=str, help='string\tinput bed file.  e.g. /in/sample.bed')
parser.add_argument(posList[2], type=str, help='string\tsamtools.  e.g. /path/to/samtools')
parser.add_argument(posList[3], type=str, help='string\toutput regional coverage statistics.  e.g. /out/sample.bed')
# optional arguments
parser.add_argument('-q', '--'+optList[0], type=int, metavar='', default=20, help='integer\tminimum mapQ for counted reads, default 20')

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


#--------------------------------------------------------
# main program
#--------------------------------------------------------
if __name__ == "__main__":

    #### coverage statistics over target regions

    bedDict = bedConfigFunc(inBed)

    out = open(outBed, 'w')
    #out.write( '\t'.join( ['chr','start','end','length','depthTotal','depthAve'] ) + '\n' ) # header for outBed

    for chr,beds in bedDict.items():

        for bed in beds:

            start,end = bed
            sys.stdout.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - INFO - ' + 'handling bed %s:%s-%s\n' %(chr,start,end)); sys.stdout.flush()

            # obtain depth stats by samtools depth
            os.system('%s depth %s -r %s:%s-%s -Q %s > %s.tmp' %(samtools,inBam,chr,start,end,mapQCut,outBed))
            depthDict = {}
            with open(outBed+'.tmp') as f:
                for line in f:
                    depthDict[ int(line.strip().split('\t')[1]) ] = int(line.strip().split('\t')[2])

            #### statistics
            for head in range(int(start),int(end),100):
                tail = head + 99
                if tail > int(end): tail = int(end)
                segTotal = sum( [ depthDict[j] for j in range(head,tail+1) if j in depthDict ] )
                segLength = tail - head + 1
                segAve = segTotal / segLength

                out.write( '\t'.join( map(str,[chr,head,tail,format(segAve,'.2f')]) ) + '\n' )
            os.system('rm %s.tmp' %(outBed))

    out.close()


sys.stdout.write('\n' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' - ' + sys.argv[0] + ' - End of Program\n'); sys.stdout.flush()

