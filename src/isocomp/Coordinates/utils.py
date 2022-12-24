"""Stand alone functions for dealing with genomic coordinates in isocomp"""

# std lib import
import logging
import os
# ext dependencies
import pyranges as pr

logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = ['create_comparison_windows']


def create_comparison_windows(gtf_list: list, **kwargs) -> pr.PyRanges:
    """Read in gtf files to pyrange objects, concat multiple gtf files 
    together, and merge overlapping ranges

    Args:
        gtf_list (list): a list of paths to gtf files
        **kwargs (dict): optional keyword arguments for pr.merge()

    Raises:
        IOError: raised when gtf_list is not a list
        FileNotFoundError: raised when items in gtf_list are not found
        AssertionError: raised when the ext of gtf_list items is not .gtf 

    Returns:
        pr.PyRanges: A pyranges object with columns Chromosome, Start, End and 
        Strand (included by default, but merge can be performed by ignoring 
        strand by passing the appropriate pranges.merge keyword argument). Note 
        that the overlap required to merge (by default, 0) may be controlled 
        with the 'slack' keyword argument to pyranges.merge().
    """
    # check input
    logging.debug(gtf_list)
    if not isinstance(gtf_list, list):
        raise IOError('pyranges_list must be type list')
    for path in gtf_list:
        if not os.path.exists(path):
            raise FileNotFoundError(f'{path} does not exist')
        elif os.path.splitext(path)[1] != '.gtf':
            raise AssertionError(f'File extension of {path} is not gtf. This '
                                 f'must be a gtf file, and the way it is '
                                 f'confirmed is by the extension. '
                                 f'Either reformat to gtf or rename.')

    # read in data
    pyranges_list = [pr.read_gtf(x) for x in gtf_list]

    # concat ranges
    concat_ranges = pr.concat(pyranges_list)

    # merge ranges. pass any additional keyword arguments
    # from the function call to pr.merge
    merged_ranges = concat_ranges.merge(**kwargs)

    logging.debug('number of merged ranges: %s', str(len(merged_ranges)))

    return merged_ranges
