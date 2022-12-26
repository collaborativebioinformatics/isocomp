import operator
import logging

import pyranges as pr

from .IsoformLibrary import IsoformLibrary

logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = ['vector_crosser', 'align_isoforms']


# TODO this isn't very efficient b/c of the list operations. There is likely a
# better implementation maybe in numpy, or a more pythonic way of doing this
# same thing in a lot less code


def vector_crosser(v1: list, v2: list, equals: bool = False) -> dict:
    """given two lists with any length and any element type, generate a  
    a dictionary with keys 'V1' and 'V2', each of which stores a list. 
    Indicies of the list correspond to one another which describe all 
    unique combinations of the elements of v1 and v2. 
    Set equals to TRUE to return corresponding elements with equal values, 
    eg 1 == 1. This is based on R code here:
    https://github.com/mhesselbarth/suppoRt/blob/HEAD/R/expand_grid_unique.R

    Args:
        v1 (list): a list of items
        v2 (list): a list of items
        equals (bool, optional): whether to return paired elements where 
        the values of v1 and v2 are the same, eg '1' '1' would be in the same 
        index in V1 and V2 if this is set to True. Defaults to False.

    Returns:
        dict: a dictionary with keys 'V1' and 'V2', each of which stores a 
        list. Indicies of the list correspond to one another which describe 
        all unique combinations of the elements of v1 and v2
    """
    d = {}

    unique_v1 = list(set(v1))
    unique_v2 = list(set(v2))

    def inner(i: int) -> None:
        """This is intended to be used in the for loop below. The variable 
        z stores the set diff between unique_v2 and, depending on the value 
        of i and the variable equals, some range of unique_v1. For example, 
        in the for loop below, we iterate over the length of unique_v1. If the 
        length is three, __and__ equals is set to False, then the first 
        iteration takes the set diff of unique_v2 and unique_v1[0:1] which 
        is the first element of unique_v1. If equals is set to True, then the 
        first iteration is the set diff of unique_v2 and unique_v1[0:0] which 
        returns the entirety of unique_v2. this continues in the for loop below, 
        iteratively taking more of unique_v1

        Args:
            i (int): This is used to extract a range of unique_v1 in the 
            set difference operation, and to extract a a given value from 
            unique_v1 and append it (repeated for length(z)) to V1 while 
            z (the set diff result) is appended to V2
        """
        z = list(set(unique_v2) - set(unique_v1[0:i + operator.not_(equals)]))
        if z:
            d.setdefault('V1', []).extend([unique_v1[i]]*len(z))
            d.setdefault('V2', []).extend(z)

    # see the docstring for inner() above
    for i in range(len(unique_v1)):
        inner(i)

    return d


def align_isoforms(isoform1, isoform2) -> dict:

    aln = edlib.align(
        isoform1,
        isoform2,
        mode='NW',
        task='path')

    # TODO figure out this 'normalization' -- should it be the longest
    # iso in the cluster? longest between iso1 and iso2? Does it matter?
    out = {'normalized_edit_dist':
           round(aln['editDistance'] / len(isoform1), 2),
           'cigar': aln['cigar']}

    return out

def compare_isoforms_in_cluster(
    isoform_library: IsoformLibrary, cluster: int) -> dict:
    
    if cluster not in isoform_library.cluster_list:
        raise ValueError(f'{cluster} not in isoform_library.cluster_list')

    cluster_gtf = isoform_library.clustered_gtf


# def find_unique_isoforms(
#     transcripts: pysam.AlignmentFile, clustered_gtf: pr.PyRanges):

# def inner(gtf_cluster: pr.PyRanges):

#     align_res_dict = {}
#     normalized_edit_dist_list = []

    # add unique identifier to gtf

    # create cross

    # iterate over cross, extracting corresponding seqs and
    # send running align_isoform(). Add result to align_res_dict

    # select alignments by percentile edit distance
    # if alns:
    #     cut = numpy.percentile(eds, 100 - minPercent)
    #     picks = ['__'.join([str(ed), allNames[i], allContigs[i], str(
    #         allStarts[i]), allSeqs[i], cigars[i]]) for i, ed in enumerate(eds) if ed <= cut]
    # else:
    #     picks = []

# # function to config input bed & bam & fa files
# def readSeqsFunc(inBed, manifest):

#     '''
#     read all input information.

#     @parameter inBed - input bed file in "chr start end ..." format.
#     @parameter manifest - input manifest file in "sample_id,path_to_bam,path_to_fasta" format.

#     '''

#     outDict = {}
#     '''
#     {
#     locus1: {
#         sample1: {
#             qname1: {
#                 seq:
#                 start:
#                 length:
#             }
#             qname2: {...}
#         }
#         sample2: {...}
#     locus2: {...}
#     ...
#     }
#     '''

#     # read windows
#     with open(inBed) as f:
#         for l in f: outDict[tuple(l.strip().split()[:3])] = {}

#     # read manifest
#     with open(manifest) as f: fields = [ l.strip().split(',') for l in f ]

#     # read sequences
#     for sample,bam,fa in fields:

#         # read fa & bam
#         seqDict = {s.id : str(s.seq) for s in SeqIO.parse(fa, 'fasta')}
#         samfile = pysam.AlignmentFile(bam, 'rb')

#         # config sequences
#         for locus in outDict.keys():

#             outDict[locus][sample] = {}

#             for read in samfile.fetch(locus[0], int(locus[1]), int(locus[2])):
#                 qname = read.query_name.split('_')[2]
#                 outDict[locus][sample][qname] = {}
#                 outDict[locus][sample][qname]['seq'] = seqDict[qname]
#                 outDict[locus][sample][qname]['contig'] = locus[0]
#                 outDict[locus][sample][qname]['start'] = read.reference_start + 1 # convert 0-based to 1-based
#                 outDict[locus][sample][qname]['query_length'] = read.query_length

#         samfile.close()

#     return outDict


# # function for unique isoform search and N-W ed stats
# def compareFunc(seqDict, outFile, minPercent):

#     '''
#     find unique isoforms and corresponding N-W alignment stats.

#     @parameter seqDict - configed sequences from readSeqsFunc().
#     @parameter outFile - output file for results.
#     @parameter minPercent - minimum percentile for edit distance output (descending order).

#     '''

#     out = open(outFile, 'w')
#     out.write('\t'.join(['win_chr',
#                          'win_start',
#                          'win_end',
#                          'total_isoform',
#                          'isoform_name',
#                          'sample_from',
#                          'sample_compared_to',
#                          'mapped_start',
#                          'isoform_sequence',
#                          'selected_alignments']) + '\n')

#     for locus,allSamples in seqDict.items():

#         # all isoforms of the locus
#         '''allSeqs, allNames = [], []
#         for s,d in allSamples.items():
#             allSeqs += [ d[q]['seq'] for q in d.keys() ]
#             allNames += [ s+'_'+q for q in d.keys() ]
#         total = len(allSeqs)
#         '''
#         total = sum([ len([ allSamples[sample][qname]['seq'] for qname in allSamples[sample].keys() ]) for sample in allSamples.keys() ])

#         lookup = {}

#         for sample,sampleDict in allSamples.items():

#             # isoforms in this sample
#             this = [ sampleDict[q]['seq'] for q in sampleDict.keys() ]

#             for anotherSample in allSamples.keys():

#                 if anotherSample == sample: continue

#                 # isoforms in another sample
#                 another = [ allSamples[anotherSample][q]['seq'] for q in allSamples[anotherSample].keys() ]
#                 # intersect to obtain unique isoforms in this sample
#                 unique = set(this).difference(set(another))

#                 # N-W alignment stats
#                 for u in unique:

#                     # all condidate sequences for N-W alignment
#                     allSeqs, allContigs, allStarts, allNames = [], [], [], []
#                     for s,d in allSamples.items():
#                         allSeqs += [ d[q]['seq'] for q in d.keys() if d[q]['seq'] != u]
#                         allContigs += [ d[q]['contig'] for q in d.keys() if d[q]['seq'] != u]
#                         allStarts += [ d[q]['start'] for q in d.keys() if d[q]['seq'] != u]
#                         allNames += [ s+'_'+q for q in d.keys() if d[q]['seq'] != u ]

#                     # N-W alignment and save for quick lookup
#                     for a in allSeqs:
#                         if (u,a) not in lookup: lookup[(u,a)] = edlib.align(u, a, mode = 'NW', task = 'path')

#                     # alignment stats
#                     alns = [ lookup[(u,a)] for a in allSeqs ]
#                     eds = [ round(aln['editDistance']/len(u),2) for aln in alns ]
#                     cigars = [ aln['cigar'] for aln in alns ]

#                     # select alignments by percentile edit distance
#                     if alns:
#                         cut = numpy.percentile(eds, 100 - minPercent)
#                         picks = [ '__'.join([str(ed),allNames[i],allContigs[i],str(allStarts[i]),allSeqs[i],cigars[i]]) for i,ed in enumerate(eds) if ed <= cut ]
#                     else:
#                         picks = []

#                     # output
#                     for qname,v in sampleDict.items():
#                         temp = [total, qname, sample, anotherSample, v['start'], u, ','.join(picks)]
#                         if v['seq'] == u: out.write('\t'.join( map(str, list(locus)+temp) ) +'\n')

#     out.close()
