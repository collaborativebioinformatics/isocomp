import logging
import os
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import pandas as pd
from .IsoformLibrary import IsoformLibrary
from .compare_isoforms_in_cluster import compare_isoforms_in_cluster

# uncomment this when output filtering is implemented
#from .filter_comparisons import filter_comparisons

logger = logging.getLogger(__name__)

__all__ = ['find_unique_isoforms']


def process_cluster(cluster: str,
                    clustered_gtf: str,
                    fasta_dict: dict) -> dict:
    """Process a cluster in parallel.

    Args:
        cluster (str): The cluster ID to process
        clustered_gtf (str): path to clustered_regions.gtf
        fasta_dict (dict): A dictionary where the key is one of the
        factor levels of the cluster_regions.gtf Source column and the value
        is a path to a fasta file which stores the isoform sequences

    Returns:
        dict: A dictionary containing detailed information about the
            comparison of the two isoforms, including the cluster ID,
            chromosome, information about each isoform, and alignment details.
    """
    # Create IsoformLibrary within each process
    il = IsoformLibrary(clustered_gtf, fasta_dict)
    cluster = str(cluster)
    return compare_isoforms_in_cluster(il, cluster)


def find_unique_isoforms(clustered_gtf: str,
                         fasta_dict: dict,
                         num_cores=None) -> pd.DataFrame:
    """Iterate over the clusters in clustered_gtf. Compare isoforms
    within clusters.

    Args:
        clustered_gtf (str): path to clustered_regions.gtf
        fasta_dict (dict): A dictionary where the key is one of the
        factor levels of the cluster_regions.gtf Source column and the value
        is a path to a fasta file which stores the isoform sequences
        num_cores (int): The number of cores to use for parallel processing.

    Returns:
        DataFrame: A dataframe which describes the isoforms which are less than
        the min_percentile similar to the other isoforms in its bin
    """
    all_comparisons = []

    # Check available CPUs
    available_cpus = os.cpu_count()

    # Validate max_workers
    if num_cores is None or num_cores > available_cpus:
        max_workers = max(1, available_cpus - 1)
    else:
        max_workers = num_cores

    il = IsoformLibrary(clustered_gtf, fasta_dict)

    # Use 'partial' to create a new function with necessary parameters
    func = partial(process_cluster, clustered_gtf=clustered_gtf,
                   fasta_dict=fasta_dict)

    # Parallel processing of clusters
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(func, il.cluster_list))

    # Flatten the list of lists to a single list
    for sublist in results:
        all_comparisons.extend(sublist)

    # note -- implement user input filtering here on what to return
    # this was my previous (buggy) implementation:
    # compare_df_fltr = filter_comparisons(all_comparisons)
    # return pd.DataFrame(all_comparisons) #compare_df_fltr

    # return raw result
    return pd.DataFrame(all_comparisons)