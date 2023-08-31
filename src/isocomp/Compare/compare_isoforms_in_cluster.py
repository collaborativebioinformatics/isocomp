import logging

from isocomp.Coordinates import Window
from .vector_crosser import vector_crosser
from .align_isoforms import align_isoforms
from .IsoformLibrary import IsoformLibrary

logger = logging.getLogger(__name__)

__all__ = ['compare_isoforms_in_cluster']


def __output_dict(cluster_id: str,
                  chr: str,
                  isoform1_window: Window,
                  isoform2_window: Window = None,
                  aln: dict = None) -> dict:
    """
    Create a dictionary containing the information regarding unique isoforms.
    If in a given cluster there are sequence comparisons, then the output
    attributes for both isoform1 and isoform2 will be not empty strings, and
    there will be alignment info. Otherwise, there will be only isoform1
    details.

    Args:
        cluster_id (str): The cluster ID for the isoforms being compared.
        chr (str): The chromosome in which the isoforms are located.
            isoform1_window (Window): A Window object containing information
        about the first isoform.
        isoform2_window (Window): A Window object containing information
        about the second isoform.
        aln (dict): A dictionary containing alignment information between
            the two isoforms.
            It should include 'normalized_edit_dist' (float) and
            'cigar' (str) keys.

    Returns:
        dict: A dictionary containing detailed information about the
            comparison of the two isoforms, including the cluster ID,
            chromosome, information about each isoform, and alignment details.
    """
    cluster_id = str(cluster_id)
    chr = str(chr)

    compare_dict = {
        'cluster': cluster_id,
        'chr': chr
    }

    for window_attr in ['source', 'name', 'start', 'end', 'strand']:
        try:
            compare_dict.update({
                'isoform1_'+window_attr: getattr(isoform1_window,
                                                 window_attr)
            })
        except AttributeError as exc:
            raise AttributeError(f'isoform1_window must be a Window '
                                 f'object with the following attributes: '
                                 f'{window_attr}') from exc
        if isoform2_window:
            try:
                compare_dict.update({
                    'isoform2_'+window_attr: getattr(isoform2_window,
                                                     window_attr)
                })
            except AttributeError as exc:
                raise AttributeError(f'isoform2_window must be a Window '
                                     f'object with the following attributes: '
                                     f'{window_attr}') from exc

    if aln:
        compare_dict.update({
            'normalized_edit_dist': aln['normalized_edit_dist'],
            'cigar': aln['cigar']
        })

    return compare_dict


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

    # note that the score attribute stores the number of isoforms in
    # the window
    cluster_window = isoform_library.get_cluster_coord(cluster)

    # if there is only one individual in the cluster, then report all
    # isoforms in that cluster as unique
    if cluster_window.score < 2:
        for tx_id in cluster_gtf.transcript_id:
            isoform1_window = isoform_library.get_isoform_coord(tx_id)
            out.append(__output_dict(cluster,
                                     cluster_window.chr,
                                     isoform1_window))
    # else there are mutiple subjects -- do an all by all comparison of the
    # isoforms in the cluster
    # TODO parameterize the cases in which isoforms are compared --eg, 
    # same strand, overlap threshold, different subjects
    else:
        # group transcripts by coordinates; return unique
        cluster_gtf_grouped = cluster_gtf.df\
            .groupby(by=['Start', 'End', 'Strand'], as_index=True)

        for group, cluster_gtf_unique in cluster_gtf_grouped:
            if len(cluster_gtf_unique) > 1:        
                # this produces a cartesian product of sorts... looks something
                # like this:
                # vector_crosser(['tx_1','tx_2','tx_3'],['tx_1','tx_2','tx_3'])
                # {'V1': ['tx_2', 'tx_2', 'tx_1'], 'V2': ['tx_1', 'tx_3', 'tx_3']}
                # the V1 and V2 lists will be the same length, so if you iterate over
                # the length of either list and compare the elements at the same index,

                cross_isoforms = vector_crosser(
                    cluster_gtf_unique.unique_id,
                    cluster_gtf_unique.unique_id)

                # iterate over the comparisons produced by vector_crosser() and
                # conduct the sequence alignments
                for i in range(len(cross_isoforms['V1'])):

                    # get the unique_id corresponding to two comparisons in the
                    # cross_isoforms dict
                    isoform1_id = cross_isoforms['V1'][i]
                    isoform2_id = cross_isoforms['V2'][i]

                    # create window objects which describe the location of the isoforms
                    # according to the gtf
                    isoform1_window = isoform_library\
                        .get_isoform_coord(unique_id=isoform1_id)
                    isoform2_window = isoform_library\
                        .get_isoform_coord(unique_id=isoform2_id)

                    # compare the isoform sequences
                    aln = align_isoforms(
                        isoform_library.get_isoform_seq(unique_id=isoform1_id),
                        isoform_library.get_isoform_seq(unique_id=isoform2_id))

                    # append the compare_dict as an element to the list out
                    out.append(__output_dict(cluster,
                                            cluster_window.chr,
                                            isoform1_window,
                                            isoform2_window,
                                            aln))
            else:
                tx_id = cluster_gtf_unique['unique_id'].iloc[0]
                isoform1_window = isoform_library.get_isoform_coord(unique_id=tx_id)
                out.append(__output_dict(cluster,
                                            cluster_window.chr,
                                            isoform1_window))
    return out
