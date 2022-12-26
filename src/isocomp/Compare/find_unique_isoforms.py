import logging

from pandas import DataFrame

from .IsoformLibrary import IsoformLibrary
from .compare_isoforms_in_cluster import compare_isoforms_in_cluster
from .filter_comparisons import filter_comparisons

logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = ['find_unique_isoforms']


def find_unique_isoforms(clustered_gtf: str,
                         fasta_dict: dict) -> DataFrame:

    il = IsoformLibrary(clustered_gtf, fasta_dict)

    all_comparisons = []

    for cluster in il.cluster_list:
        cluster = str(cluster)
        # only compare if there are more than 1 isoforms in the window
        if il.get_cluster_coord(cluster).score > 1:
            all_comparisons\
                .extend(compare_isoforms_in_cluster(il, cluster))

    compare_df_fltr = filter_comparisons(all_comparisons)

    return compare_df_fltr
