# pylint:disable=W0401
from importlib.metadata import version
import math

import pandas as pd

from isocomp import Compare
from .conftest import *


def test_vector_crosser():

    v1 = ['tx_'+str(x) for x in range(5)]

    cross_res = Compare.vector_crosser(v1, v1)

    assert len(cross_res['V1']) == len(cross_res['V2'])
    # length of the cross should be n C 2 where n is length of input
    # note that this is true when lists of the same length are passed, which
    # is the use case in the codebase
    assert math.comb(len(v1), 2) == len(cross_res['V1'])


def test_IsoformLibrary(clustered_gtf, fasta_dict):

    source = 'hg002_sqanti_fltr'
    isoform = 'PB.17.2'

    il = Compare.IsoformLibrary(clustered_gtf, fasta_dict)

    # test retrieving isoform sequence
    actual_subseq = il.get_isoform_seq(
        source,
        isoform,
        start=0,
        end=10)

    expected_subseq = 'GAGAGGCAGC'

    assert actual_subseq == expected_subseq

    # test that the gtf isoform and fasta isoform have the same length
    isoform_window = il.get_isoform_coord(source, isoform)

    assert il.fasta_dict[source].get_reference_length(isoform) == \
        len(isoform_window)


def test_align_isoforms(clustered_gtf, fasta_dict):

    il = Compare.IsoformLibrary(clustered_gtf, fasta_dict)

    # both from cluster 1
    isoform1 = il.get_isoform_seq('hg004_sqanti_fltr', 'PB.13.2')
    isoform2 = il.get_isoform_seq('hg005_sqanti_fltr', 'PB.17.1')

    actual = Compare.align_isoforms(isoform1, isoform2)

    assert isinstance(actual, dict)
    assert isinstance(actual['normalized_edit_dist'], float)
    assert isinstance(actual['cigar'], str)


def test_compare_isoforms_in_cluster(clustered_gtf, fasta_dict):
    il = Compare.IsoformLibrary(clustered_gtf, fasta_dict)

    cluster_window = il.get_cluster_coord(str(1))
    cluster_compare = Compare.compare_isoforms_in_cluster(il, str(1))

    # this should be the same length as the crossed vectors, which is the
    # number of tx in the window choose 2. the cluster_window.score attr
    # stores the number of tx in the window
    assert math.comb(cluster_window.score, 2) == len(cluster_compare)


def test_filter_comparisons(clustered_gtf, fasta_dict):

    # note that this code is the same as in find_unique_isoforms, but 
    # is repeated here to get all_comparisons for the asserts below

    il = Compare.IsoformLibrary(clustered_gtf, fasta_dict)
    all_comparisons = []
    for cluster in il.cluster_list:
        cluster = str(cluster)
        # only compare if there are more than 1 isoforms in the window
        if il.get_cluster_coord(cluster).score > 1:
            all_comparisons\
                .extend(Compare.compare_isoforms_in_cluster(il, cluster))

    compare_df_fltr = Compare.find_unique_isoforms(clustered_gtf, fasta_dict)

    assert len(compare_df_fltr) > 0
    assert len(compare_df_fltr) < len(pd.DataFrame(all_comparisons))
