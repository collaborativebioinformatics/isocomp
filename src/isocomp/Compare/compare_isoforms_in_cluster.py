import logging

from .vector_crosser import vector_crosser
from .align_isoforms import align_isoforms
from .IsoformLibrary import IsoformLibrary

logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = ['compare_isoforms_in_cluster']


def compare_isoforms_in_cluster(
        isoform_library: IsoformLibrary, cluster: str) -> list:
    """Given an IsoformLibrary object and a cluster identifier, compare 
    all isoforms in the cluster against one another.

    Args:
        isoform_library (IsoformLibrary): An IsoformLibrary object
        cluster (str): A cluster identifier

    Raises:
        ValueError: Raised if a given cluster has less than 2 unique isoforms

    Returns:
        list: A list of dictionaries with the same keys which may be used 
              to create a dataframe describing the isoform comparisons
    """

    # TODO this isn't the right comparison -- there are a few questions, like
    # is the comparison x1 to x2 different than x2 to x1? and what normalized
    # means. But, this should be enough to show how to use the objects to
    # easily extract the information and perform the comparisons within cluster
    
    # this will be returned from the function. The intention is that it will 
    # be filled with dictionary elements. See compare_dict in the for loop 
    # below
    out = []

    # subset the gtf down to a single cluster
    # TODO consider doing this as a view of the isoform_library
    cluster_gtf = isoform_library.get_cluster(cluster)
    
    # note that the score attribute stores the number of isoforms in the window
    cluster_window = isoform_library.get_cluster_coord(cluster)
    
    # raise error if there are fewer than 2 isoforms in the cluster
    if cluster_window.score < 2:
        raise ValueError('cluster %s has less than 2 transcripts' % cluster) #pylint:disable=C0209 # noqa:E501,E261,E262

    # this produces a cartesian product of sorts... looks something like this:
    # vector_crosser(['tx_1','tx_2','tx_3'],['tx_1','tx_2','tx_3'])
    # {'V1': ['tx_2', 'tx_2', 'tx_1'], 'V2': ['tx_1', 'tx_3', 'tx_3']}
    # see vector_crosser() docs
    cross_isoforms = vector_crosser(
        cluster_gtf.unique_id,
        cluster_gtf.unique_id)
    
    # iterate over the comparisons produced by vector_crosser() and conduct 
    # the sequence alignments
    for i in range(len(cross_isoforms['V1'])):

        # get the unique_id corresponding to two comparisons in the 
        # cross_isoforms dict
        isoform1_id = cross_isoforms['V1'][i]
        isoform2_id = cross_isoforms['V2'][i]
        
        # create window ojects which describe the location of the isoforms
        # according to the gtf
        isoform1_window = isoform_library\
            .get_isoform_coord(unique_id=isoform1_id)
        isoform2_window = isoform_library\
            .get_isoform_coord(unique_id=isoform2_id)
        
        # compare the isoform sequences
        aln = align_isoforms(
            isoform_library.get_isoform_seq(unique_id=isoform1_id),
            isoform_library.get_isoform_seq(unique_id=isoform2_id)
        )
        
        # output the result of the alignment
        compare_dict = {
            'cluster': cluster,
            'chr': cluster_window.chr,

            'isoform1_source': isoform1_window.source,
            'isoform1_name': isoform1_window.name,
            'isoform1_start': isoform1_window.start,
            'isoform1_end': isoform1_window.end,
            'isoform1_strand': isoform1_window.strand,

            'isoform2_source': isoform2_window.source,
            'isoform2_name': isoform2_window.name,
            'isoform2_start': isoform2_window.start,
            'isoform2_end': isoform2_window.end,
            'isoform2_strand': isoform2_window.strand,

            'normalized_edit_dist': aln['normalized_edit_dist'],
            'cigar': aln['cigar']
        }
        
        # append the compare_dict as an element to the list out
        out.append(compare_dict)

    return out
