# std lib import
import logging
import os
# external dependencies
import pyranges as pr
import pandas as pd
import pysam

logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = ['IsoformLibrary']


class IsoformLibrary():

    _clustered_gtf_path = None
    _clustered_gtf = None
    _cluster_list = None
    _fasta_dict = None

    def __init__(self,
                 clustered_gtf_path: str,
                 fasta_dict: dict,
                 __validate_strict: bool = True) -> None:
        # note that there is a good deal of validation of paths, data types,
        # that the Sources identified in the cluster_gtf match the fasta
        # sources (keys of the fasta_dict), and that the gtf transcripts
        # are all in the corresponding fasta
        self.clustered_gtf_path = clustered_gtf_path
        self.fasta_dict = fasta_dict
        self.__validate(__validate_strict)

    # properties getters/setters ----------------------------------------------

    @property
    def clustered_gtf_path(self) -> str:
        return self._clustered_gtf_path

    @clustered_gtf_path.setter
    def clustered_gtf_path(self, new_path: str) -> None:
        logging.debug('trying to set new clustered_gtf_path: {new_path}')
        if not os.path.exists(new_path):
            raise FileNotFoundError(f'{new_path} does not exist')
        # TODO allow gff and check format, not extension
        if os.path.splitext(new_path)[1] != ".gtf":
            raise IOError('the gtf path must end with ".gtf"')

        # open the new gtf
        self.clustered_gtf = self.__open_gtf(new_path)
        self.cluster_list = list(set(self.clustered_gtf.Cluster))

    @property
    def clustered_gtf(self) -> pr.PyRanges:
        return self._clustered_gtf

    @clustered_gtf.setter
    def clustered_gtf(self, new_gtf_pr: pr.PyRanges) -> None:
        if not isinstance(new_gtf_pr, pr.PyRanges):
            raise ValueError('clustered_gtf must be a '
                             'pyranges.PyRanges object')
        if 'unique_id' in new_gtf_pr.columns:
            if not len(set(new_gtf_pr.unique_id)) == len(new_gtf_pr):
                raise NameError('There exists a column called unique_id '
                                'but it is not unique. Either fix this or '
                                'drop the column and try again')
            self._clustered_gtf = new_gtf_pr
        else:
            self._clustered_gtf = new_gtf_pr\
                .insert(
                    pd.DataFrame.from_dict(
                        {'unique_id': ['tx_'+str(x)
                                       for x in range(len(new_gtf_pr))]}))

    @property
    def cluster_list(self) -> list:
        return self._cluster_list

    @cluster_list.setter
    def cluster_list(self, new_list: list) -> None:
        self._cluster_list = new_list

    @property
    def fasta_dict(self) -> dict:
        return self._fasta_dict

    @fasta_dict.setter
    def fasta_dict(self, new_fasta_dict: dict) -> None:

        logging.debug("new dict: %s", new_fasta_dict)

        # check type
        if not isinstance(new_fasta_dict, dict):
            raise ValueError('fasta_dict must be a dictionary where keys '
                             'correspond to the source column of the '
                             'clustered_gtf and values are paths to fasta '
                             'of isoform sequences')

        # make sure  paths exists
        for k, v in new_fasta_dict.items():
            if not os.path.exists(v):
                raise FileNotFoundError(f'file {v} for sample {k} '
                                        f'does not exist')
            if not os.path.exists(v+'.fai'):
                raise FileNotFoundError(f'index file for {v} does not exist. '
                                        'run samtools faidx to produce an '
                                        'index file and try again.')

        # update fasta_dict
        self._fasta_dict = self.__open_fastas(new_fasta_dict)

    # private methods ---------------------------------------------------------

    def __open_gtf(self, gtf_path: str) -> pr.PyRanges:
        return pr.read_gtf(gtf_path)

    def __open_fastas(self, fasta_dict) -> dict:
        return {k: pysam.FastaFile(v) for k, v in fasta_dict.items()}  # pylint:disable=E1101 # noqa: E501

    def __validate(self, strict: bool = True) -> None:

        # validate that the set of 'Sources' in the source column and the
        # keys of fasta_dict are the same
        set_diff_source_1 = \
            set(self.clustered_gtf.Source) - set(self.fasta_dict.keys())
        set_diff_source_2 = \
            set(self.fasta_dict.keys()) - set(self.clustered_gtf.Source)

        source_msg = 'The following values are in the %s but not in %s: %s'

        if set_diff_source_1:
            msg_tmp = source_msg % (
                'gtf Source column', 'fasta dict keys', set_diff_source_1)
            if strict:
                raise ValueError(msg_tmp)
            else:
                Warning(msg_tmp)

        elif set_diff_source_2:
            msg_tmp = source_msg % (
                'fasta dict keys', 'gtf Source column', set_diff_source_1)
            if strict:
                raise ValueError(msg_tmp)
            else:
                Warning(msg_tmp)

        # validate that the fasta file corresponding to a given source
        # has all of the transcripts in it
        isoform_msg = 'the following values are in the %s gtf but not in ' +\
            'the corresponding fasta file: %s'
        for src in self.fasta_dict.keys():
            # TODO handle hardcoding on transcript_id -- this should be
            # configurable
            tx_list = set(self.clustered_gtf[self.clustered_gtf.Source == src]
                          .transcript_id)
            set_diff_isoforms = tx_list - set(self.fasta_dict[src].references)

            if set_diff_isoforms:
                msg_tmp = isoform_msg % (src, set_diff_isoforms)
                if strict:
                    raise ValueError(msg_tmp)
                else:
                    Warning(msg_tmp)

    # public methods ----------------------------------------------------------

    def get_isoform(self, source: str, isoform: str, **kwargs) -> str:

        if source not in self.fasta_dict.keys():
            raise KeyError(f'the source: {source} is not in the '
                           f'fasta_dict keys: {self.fasta_dict.keys()}')
        if isoform not in set(self.clustered_gtf.transcript_id):
            raise ValueError(f'the isoform {isoform} is not in the '
                             f'transcript_id columns of clustered_gtf')

        return self.fasta_dict[source].fetch(reference=isoform, **kwargs)

    def get_cluster(self, cluster: int) -> pr.PyRanges:

        if cluster not in self.cluster_list:
            raise ValueError(f'{cluster} not in cluster_list')

        # TODO remove hard coded Cluster to external config
        return self.clustered_gtf[self.clustered_gtf.Cluster == cluster]
