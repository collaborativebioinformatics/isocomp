import logging

import pandas as pd

logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = ['filter_comparisons']


def filter_comparisons(comparison_list: list,
                       min_quantile: float = 0.05) -> pd.DataFrame:

    df = pd.DataFrame(comparison_list)

    cut = df.normalized_edit_dist.quantile(min_quantile)

    df_fltr = df[df.normalized_edit_dist <= cut]

    return df_fltr
