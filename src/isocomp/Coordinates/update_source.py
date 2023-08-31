# pylint:disable=W0108
# std lib import
import logging
# ext dependencies
import pyranges as pr

logger = logging.getLogger(__name__)

__all__ = ['update_source']


def update_source(gtf_pr: pr.PyRanges, new_source: str) -> pr.PyRanges:
    """Update the value in the source column of a pyranges obj read in 
    from a gtf/gff

    Args:
        gtf_pr (pr.PyRanges): a py.PyRanges object read in from a gff/gtf
        new_source (str): the value to enter into the Source column

    Returns:
        pr.PyRanges: the input pyranges obj with the updated Source column
    """
    if 'Source' not in gtf_pr.columns:
        raise KeyError(f'Column "Source" does not exist in '
                       f'{new_source+".gtf/gff"}')

    # update the value in the Source column
    gtf_pr.Source = new_source

    return gtf_pr