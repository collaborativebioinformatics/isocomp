import logging
import pandas as pd
from pandas import DataFrame

from .IsoformLibrary import IsoformLibrary
from .compare_isoforms_in_cluster import compare_isoforms_in_cluster
from .filter_comparisons import filter_comparisons

logger = logging.getLogger(__name__)

__all__ = ['find_unique_isoforms']


def find_unique_isoforms(clustered_gtf: str,
                         fasta_dict: dict) -> DataFrame:
    """Iterate over the clusters in clustered_gtf. Compare isoforms 
    within clusters.

    Args:
        clustered_gtf (str): path to clustered_regions.gtf
        fasta_dict (dict): A dictionary where the key is one of the 
        factor levels of the cluster_regions.gtf Source column and the value 
        is a path to a fasta file which stores the isoform sequences

    Returns:
        DataFrame: A dataframe which describes the isoforms which are less than 
        the min_percentile similar to the other isoforms in its bin
    """
    # instantiate an IsoformLibrary
    il = IsoformLibrary(clustered_gtf, fasta_dict)
    # instantiate a list to store the dict objects which result from 
    # running the iteration below
    all_comparisons = []
    # iterate over clusters and compare isoforms
    for cluster in il.cluster_list:
        cluster = str(cluster)
        logger.debug(cluster)
        # only compare if there are more than 1 isoforms in the window
        if il.get_cluster_coord(cluster).score > 1:
            all_comparisons\
                .extend(compare_isoforms_in_cluster(il, cluster))
    # filter the result of the comparisons
    #compare_df_fltr = filter_comparisons(all_comparisons)

    return pd.DataFrame(all_comparisons) #compare_df_fltr
