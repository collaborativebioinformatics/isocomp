# std lib
import sys
import os
import logging
from logging.config import dictConfig
import argparse
from typing import Callable
from importlib.metadata import version
# local imports
from .Coordinates import utils

logging.getLogger(__name__).addHandler(logging.NullHandler())


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
        'create_windows': "from a list of SQANTI gtfs, create a bed6 "
        "formatted file -- 0 indexed, half open intervals -- where each entry "
        "in the bed file represents a discrete region of coverage over the "
        "reference genome."
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

    # parse_bam subparser -----------------------------------------------------
    subparsers = parser.add_subparsers(
        help="Available Tools")

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

    create_windows_output = create_windows_parser.add_argument_group('output')
    create_windows_output.add_argument(
        "-p",
        "--output_prefix",
        help="path to output directory. if not provided, " +
        "output to current working directory",
        default="",
        required=False)
    create_windows_output.add_argument(
        "--overwrite",
        help="By default, create_windows will overwrite a file with the same " +
        "name in the output directory. Set --no-overwrite to avoid overwriting.",
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

    # return the top level parser to be used in the main method below ---------
    return parser


def __create_windows_gtfs(args=None):
    """Entry point to merge gtf regions and write a merged gtf file. 
    Expected arguments are input(list of gtf paths),output_prefix(str) and 
    overwrite(bool)"""
    logging.debug(args)

    merged_regions = utils.create_comparison_windows(args.input)

    # add Name and Score columns per bed6 format
    merged_regions.Name = ['region_'+str(x)
                           for x in range(len(merged_regions))]
    merged_regions.Score = [1000 for x in range(len(merged_regions))]

    output_filename = args.output_prefix+'.bed' \
        if args.output_prefix \
        else 'merged_regions.bed'
    logging.debug(output_filename)

    if os.path.exists(output_filename) and not args.overwrite:
        raise FileExistsError(f'file with name {output_filename} already '
                              f'exists set overwrite to True if you wish '
                              'to overwrite')

    # write out as a bed6 format file
    merged_regions\
        .df[['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand']]\
        .to_csv(output_filename, sep='\t', header=None)


def main(args=None) -> None:
    """Entry point to isocomp"""

    # parse the cmd line arguments
    arg_parser = parse_args()

    args = arg_parser.parse_args(args)

    # this is a default setting -- if it is not set, it means
    # that nothing was passed on the cmd line. Instead, print the
    # help message
    try:
        log_level = args.log_level.upper()
    except AttributeError:
        sys.exit(arg_parser.print_help())

    # set the logging details
    log_config = {
        "version": 1,
        "root": {
            "handlers": ["console"],
            "level": f"{log_level}"
        },
        "handlers": {
            "console": {
                "formatter": "std_out",
                "class": "logging.StreamHandler"
            }
        },
        "formatters": {
            "std_out": {
                "format": "%(asctime)s : %(module)s : " +
                "%(funcName)s : line: %(lineno)d\n" +
                "\tprocess details : %(process)d, %(processName)s\n" +
                "\tthread details : %(thread)d, %(threadName)s\n" +
                "\t%(levelname)s : %(message)s",
                "datefmt": "%Y-%m-%d %H:%M:%S"
            }
        }
    }
    dictConfig(log_config)
    # log the cmd line arguments at the debug level
    logging.debug(sys.argv)
    logging.debug(str(args))

    # note that this works b/c the subparser set_defaults function attribute
    # is set.
    # see https://docs.python.org/3/library/argparse.html#parser-defaults
    # scroll up from that point to see a usage example
    args.func(args)


if __name__ == "__main__":
    sys.exit(main())
