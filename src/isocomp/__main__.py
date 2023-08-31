# std lib
import sys
import os
import logging
from logging.config import dictConfig
import argparse
from typing import Callable
from importlib.metadata import version
# extenal dependencies
import pandas as pd
from Bio import SeqIO
# local imports
from .Coordinates import create_comparison_windows
from .Compare import find_unique_isoforms
from .utils.configure_logging import configure_logging

logger = logging.getLogger(__name__)


def parse_args() -> Callable[[list], argparse.Namespace]:
    """Create a cmd line argument parser for isocomp
    Returns:
        Callable[[list],Namespace]: This function returns the main argparser. 
        If the subparser set_defaults method is set, it can make the correct 
        decision as to which function to call to execute a given tool. See 
        the main method below for usage
    """
    # create the top-level parser

    script_descriptions = {
        'fasta_to_fastq': "convert a fasta file to a fastq file with max "
        "quality for each base pair",
        'create_windows': "from a list of SQANTI gtfs, create a bed6 "
        "formatted file -- 0 indexed, half open intervals -- where each entry "
        "in the bed file represents a discrete region of coverage over the "
        "reference genome.",
        'find_unique_isoforms': "Using the clustered gtf output of "
        "create_windows and a csv file with columns 'source' and 'fasta' "
        "which provide a map from the unique values in 'Source' in the "
        "clustered gtf to the corresponding fasta file, find unique "
        "isoforms in the clusters described in the clustered_gtf."
    }

    # common options -- these can be applied to all scripts via the 'parent'---
    # argument -- see subparsers.add_parser for parse_bam below ---------------
    common_args = argparse.ArgumentParser(prog="isocomp", add_help=False)
    common_args_group = common_args.add_argument_group('general')
    common_args_group.add_argument(
        "-l",
        "--log_level",
        choices=("critical", "error", "warning", "info", "debug"),
        default="warning")

    # Create a top level parser -----------------------------------------------
    parser = argparse.ArgumentParser(
        prog='isocomp',
        description=f"isocomp: {version('isocomp')}")
    # add argument to get version
    parser.add_argument(
        "-v",
        "--version",
        action='version',
        version='%(prog)s '+f'{version("isocomp")}')
    
    # create a subparser
    subparsers = parser.add_subparsers(
        help="Available Tools")

    # fasta_to_fastq subparser ------------------------------------------------

    fasta_to_fastq_parser = subparsers.add_parser(
        'fasta_to_fastq',
        help=script_descriptions['fasta_to_fastq'],
        description=script_descriptions['fasta_to_fastq'],
        prog='fasta_to_fastq',
        parents=[common_args])

    # set the function to call when this subparser is used
    fasta_to_fastq_parser.set_defaults(func=__fasta_to_fastq)

    fasta_to_fastq_input = \
        fasta_to_fastq_parser.add_argument_group('input')

    fasta_to_fastq_input.add_argument(
        "-i",
        "--fasta",
        help="path to a fasta file. Output is to standard out",
        required=True)

    # create_windows subparser ------------------------------------------------
    
    create_windows_parser = subparsers.add_parser(
        'create_windows',
        help=script_descriptions['create_windows'],
        description=script_descriptions['create_windows'],
        prog='create_windows',
        parents=[common_args])

    # set the function to call when this subparser is used
    create_windows_parser.set_defaults(func=__create_windows_gtfs)

    create_windows_input = create_windows_parser.add_argument_group('input')
    create_windows_input.add_argument(
        "-i",
        "--input",
        nargs='+',
        help="a space delimited list of paths to gtf files from a " +
        "concatenated merged set of regions will be created",
        required=True)
    create_windows_input.add_argument(
        "-f",
        "--feature",
        help="The feature on which to cluster. Default is transcript",
        default="transcript",
        required=True)

    create_windows_output = create_windows_parser.add_argument_group('output')
    create_windows_output.add_argument(
        "-o",
        "--output_prefix",
        help="This should be either a file basename, which "
        "will be written to the current working directory, or a valid path "
        "up to a file basename, eg /output/path/concat would result in the "
        "file output/path/concat.gtf being written",
        default="",
        required=False)
    create_windows_output.add_argument(
        "--overwrite",
        help="By default, create_windows will overwrite a file with the same "
        "name in the output directory. "
        "Set --no-overwrite to avoid overwriting.",
        action=argparse.BooleanOptionalAction,
        required=False)

    # TODO consider adding options to adjust how pyranges.merge behaves here,
    # eg stranded and slack
    # create_windows_settings = create_windows_parser\
    #     .add_argument_group('settings')
    # create_windows_settings.add_argument(
    #     "-q",
    #     "--mapq_threshold",
    #     help="Reads less than or equal to mapq_threshold " +
    #     "will be marked as failed",
    #     type=int,
    #     default=10)

    # find_unique_isoforms subparser ------------------------------------------

    find_unique_isoforms_parser = subparsers.add_parser(
        'find_unique_isoforms',
        help=script_descriptions['find_unique_isoforms'],
        description=script_descriptions['find_unique_isoforms'],
        prog='find_unique_isoforms',
        parents=[common_args])

    # set the function to call when this subparser is used
    find_unique_isoforms_parser.set_defaults(func=__find_unique_isoforms)

    find_unique_isoforms_input = \
        find_unique_isoforms_parser.add_argument_group('input')

    find_unique_isoforms_input.add_argument(
        "-a",
        "--clustered_gtf",
        help="path to the clustered gtf file, likely created by create_windows",
        required=True)

    find_unique_isoforms_input.add_argument(
        "-f",
        "--fasta_map",
        help="path to a csv file with columns 'source' and 'fasta' which "
        "provide a map between the unique factor levels of the Source "
        "column of the clustered gtf and a corresponding fasta file which "
        "stores the transcript sequences",
        required=True)

    find_unique_isoforms_output = \
        find_unique_isoforms_parser.add_argument_group('output')
    find_unique_isoforms_output.add_argument(
        "-p",
        "--output_prefix",
        help="This should be either a file basename -- NO extension -- which "
        "will be written to the current working directory, or a valid path "
        "up to a filename, eg /output/path/concat would result in the file "
        "output/path/concat.gtf being written",
        default="",
        required=False)
    find_unique_isoforms_output.add_argument(
        "--overwrite",
        help="By default, find_unique_isoforms will overwrite a file "
        "with the same name in the output directory. "
        "Set --no-overwrite to avoid overwriting.",
        action=argparse.BooleanOptionalAction,
        required=False)

    # return the top level parser to be used in the main method below ---------
    return parser


# TODO input to this function should be the gtfs and fastas which will be used
# to do the comparisons. This should call a function which takes the seqnames
# from each of the fastas and finds the corresponding column in a pyrange
# object created from the gtf. If the gtf does not fully describe the fasta,
# throw and error. Output something a json config file suitable to configure
# the isocomp processes (eg, maybe 'transcript' isn't the feature to use from
# col3 of the gtf or some such)
def __validate_input(args=None) -> None:
    """_summary_

    Args:
        args (_type_, optional): _description_. Defaults to None.
    """
    raise NotImplementedError()


def __fasta_to_fastq(args=None) -> None:
    """Convert a fasta file to a fastq file. This function simply iterates 
    over a fasta file and prints a fastq file line for each line of the fasta 
    file

    Args:
        args (argparse.Namespace, optional): Expects one argument, named 
        fasta, a path to a fasta format file. Defaults to None.
    """

    logger.debug(args)

    for seq in SeqIO.parse(args.fasta, "fasta"):
        seq.letter_annotations["solexa_quality"] = [40] * len(seq)
        print(seq.format("fastq"), end='')


def __create_windows_gtfs(args=None) -> None:
    """Entry point to merge gtf regions and write a merged gtf file. 
    Expected arguments are input(list of gtf paths),output_prefix(str) and 
    overwrite(bool)

    Args:
        args (argparse.Namespace, optional): argparse.Namespace parsed from 
        the cmd line

    Raises:
        FileExistsError: raised if the output path exists and overwrite is 
        False
    """
    logger.debug(args)

    # TODO consider stripping extension, if one is passed, from output_prefix
    output_filename = args.output_prefix+'.gtf' \
        if args.output_prefix \
        else 'clustered_regions.gtf'
    logger.debug(output_filename)

    if os.path.exists(output_filename) and not args.overwrite:
        raise FileExistsError(f'file with name {output_filename} already '
                              f'exists set overwrite to True if you wish '
                              'to overwrite')

    clustered_regions = create_comparison_windows(args.input, args.feature)

    # write out as a bed6 format file
    clustered_regions.to_gtf(output_filename)


def __find_unique_isoforms(args=None) -> None:
    """Entry point script to compare transcripts within each cluster in the 
    clustered_gtf

    Args:
        args (argparse.Namespace, optional): argparse.Namespace parsed from 
        the cmd line

    Raises:
        FileNotFoundError: raised if the input is not found
        FileExistsError: raised if the output path exists and overwrite is 
        False
        KeyError: raised if the fasta_map doesn't conform to expectations on 
        column names
    """

    logger.debug(args)

    for path in [args.clustered_gtf, args.fasta_map]:
        if not os.path.exists(path):
            raise FileNotFoundError('%s does not exist' % path)  # pylint:disable=C0209 # noqa:E501,E261,E262

    # TODO consider stripping extension, if one is passed, from output_prefix
    output_filename = args.output_prefix+'.csv' \
        if args.output_prefix \
        else 'unique_isoforms.csv'
    logger.debug(output_filename)

    if os.path.exists(output_filename) and not args.overwrite:
        raise FileExistsError(f'file with name {output_filename} already '
                              f'exists set overwrite to True if you wish '
                              'to overwrite')

    # read in the fasta_map csv to a DataFrame and do some checking
    fasta_df = pd.read_csv(args.fasta_map)
    if 'source' not in fasta_df.columns or 'fasta' not in fasta_df.columns:
        raise KeyError('the fasta map must be a csv with labelled columns '
                       'source, which stores the sample name corresponding to '
                       'the clustered gtf Source column, and fasta, which '
                       'stores the path to a fasta file which describes the '
                       'isoform sequences')

    # create the fasta_dict from the DataFrame
    fasta_dict = dict(zip(fasta_df.source, fasta_df.fasta))

    # compare within each cluster and filter the results
    comparison_fltr_df = find_unique_isoforms(args.clustered_gtf, fasta_dict)

    # write out the results
    comparison_fltr_df.to_csv(output_filename, index=False)


def main(args=None) -> None:
    """Entry point to isocomp. Note that aside from messing with how the 
    logger is configured, you do not need to do anything with this function 
    to add more cmd line scripts"""

    # parse the cmd line arguments
    arg_parser = parse_args()

    args = arg_parser.parse_args(args)

    # this is a default setting -- if it is not set, it means
    # that nothing was passed on the cmd line. Instead, print the
    # help message
    try:
        log_level = args.log_level.upper()
        if log_level not in ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG']:
            raise ValueError("The logging level must be one of debug, "
                             "info, warning, error, "
                             "or critical.")
    except AttributeError:
        sys.exit(arg_parser.print_help())

    configure_logging(log_level)
    # log the cmd line arguments at the debug level
    logger.debug(sys.argv)
    logger.debug(str(args))

    # note that this works b/c the subparser set_defaults function attribute
    # is set.
    # see https://docs.python.org/3/library/argparse.html#parser-defaults
    # scroll up from that point to see a usage example
    args.func(args)


if __name__ == "__main__":
    sys.exit(main())
