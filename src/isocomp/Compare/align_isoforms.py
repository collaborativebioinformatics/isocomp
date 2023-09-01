# std lib
import logging
# external dependencies
import edlib

logger = logging.getLogger(__name__)

__all__ = ['align_isoforms']


def align_isoforms(isoform1: str, isoform2: str) -> dict:
    """Given two sequences, conduct needleman-wunsch alignment 
    with edlib

    Args:
        isoform1 (str): Any string, but the intent is that this is a 
        DNA string. Note that the 
        isoform2 (str): Any string, but the intent is that this is a 
        DNA string

    Returns:
        dict: A dictionary with keys 'normalized_edit_distance' and 'cigar'. 
              normalized_edit_dist is currently calculated by dividing the edit 
              distance of isoform1 vs isoform2 by the length of isoform1. 'cigar' 
              stores the cigar string.
    """

    aln = edlib.align(
        isoform1,
        isoform2,
        mode='NW',
        task='path')

    # TODO figure out this 'normalization' -- should it be the longest
    # iso in the cluster? longest between iso1 and iso2? Does it matter?
    # hi -- i am in the new_feature branch
    out = {'normalized_edit_dist':
           round(aln['editDistance'] / len(isoform1), 2),
           'cigar': aln['cigar']}

    return out
