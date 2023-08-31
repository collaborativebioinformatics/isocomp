"""Stand alone functions for dealing with genomic coordinates in isocomp"""
# pylint:disable=W0108
# std lib import
import logging
import os
# ext dependencies
import pyranges as pr

from .update_source import update_source

logger = logging.getLogger(__name__)

__all__ = ['create_comparison_windows']


def create_comparison_windows(gtf_list: list, 
                              feature: str = 'transcript', 
                              **kwargs) -> pr.PyRanges:
    """Read in gtf files to pyrange objects. To each gtf, replace the 'source' 
    with the base filename (no extention) of the gtf file. Filter on Feature 
    'transcript' and cluster overlapping ranges. Each cluster will be 
    sequentially numbered, so cluster 1 will comprise a discrete group of 
    transcripts with overlapping ranges, as will cluster 2, 3, ...

    Args:
        gtf_list (list): a list of paths to gtf files
        feature (str): The feature on which to cluster. Default is 'transcript'
        **kwargs (dict): optional keyword arguments for pr.merge()

    Raises:
        IOError: raised when gtf_list is not a list
        FileNotFoundError: raised when items in gtf_list are not found
        AssertionError: raised when the ext of gtf_list items is not .gtf 

    Returns:
        pr.PyRanges: A pyranges object with the columns Chromosome, Source 
        (source is replaced with the filename of the gtf file, which must be 
        unique in the set) Feature, Start, End, Score, Strand, Frame, 
        transcript_id, gene_id, Cluster
    """
    # check input
    logger.debug(gtf_list)
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
    pyranges_list = [update_source(pr.read_gtf(x),
                                   os.path.splitext(os.path.basename(x))[0])
                     for x in gtf_list]

    # concat gtfs
    concat_ranges = pr.concat(pyranges_list)

    # merge ranges. pass any additional keyword arguments
    # from the function call to pr.merge
    # TODO make sure that there is a way of validating that feature is a valid
    # feature in the gtf file. See
    # __validate_input() in __main__ for note on creating a config json
    concat_ranges = concat_ranges[concat_ranges.Feature == feature]
    clustered_ranges = concat_ranges.cluster(**kwargs)

    logger.debug('number of merged ranges: %s',
                  str(max(clustered_ranges.Cluster)))

    return clustered_ranges
