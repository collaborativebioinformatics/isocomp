import logging

import pandas as pd

logger = logging.getLogger(__name__)

__all__ = ['filter_comparisons']


def filter_comparisons(comparison_list: list,
                       min_quantile: float = 0.05) -> pd.DataFrame:
    """Filter the isoform comparison vector which are less than 
    min_quantile similar to each other.Note that 

    Args:
        comparison_list (list): This function is intended to be used in 
            find_unique_isoforms(), which creates a list of dict objects by 
            iterating over the clusters and storing the result of 
            compare_isoforms_incluster() in a list. This argument takes that 
            list of dict objects 
        min_quantile (float, optional): Only accept 
            isoforms which are less than the min_quantile of the normalized 
            edit distance similar to one another. Defaults to 0.05.

    Returns:
        pd.DataFrame: A filtered dataframe with the columns corresponding 
                      to the fields recorded in find_unique_isoforms()
    """
    # transform the comparison_list into a dataframe
    df = pd.DataFrame(comparison_list)
    # determine the cut point according to the normalized_edit_dist
    cut = df.normalized_edit_dist.quantile(min_quantile)
    # filter for those isoforms which are less similar than the cut point
    df_fltr = df[df.normalized_edit_dist <= cut]
    
    return df_fltr
