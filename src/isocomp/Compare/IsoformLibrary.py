# std lib import
import logging
import os
# external dependencies
import pyranges as pr
import pandas as pd
import pysam

# local imports
from isocomp.Coordinates import Window

logger = logging.getLogger(__name__)

__all__ = ['IsoformLibrary']


class IsoformLibrary():
    """This class offers a method of extracting data from a 
        clustered_regions gtf file and a corresponding fasta file with 
        isoform sequences

    Args:
        clustered_gtf_path (str): path to a clustered gtf file, intended to 
        be one output by create_windows()
        fasta_dict (dict): A dictionary where the keys are factor levels of 
        the Source column in the clustered_regions gtf and the value is 
        a path to a fasta which stores isoform sequences
        validate_strict (bool, optional): Set to False to warn, rather than 
        error, if the fasta_dict keys do not fully match the gtf Source 
        factor levels, and the transcripts in the gtf do not fully match 
        the sequence names in the fasta. Defaults to True.
    """

    _clustered_gtf_path = None
    _clustered_gtf = None
    _cluster_list = None
    _fasta_dict = None

    def __init__(self,
                 clustered_gtf_path: str,
                 fasta_dict: dict,
                 validate_strict: bool = True) -> None:
        # note that there is a good deal of validation of paths, data types,
        # that the Sources identified in the cluster_gtf match the fasta
        # sources (keys of the fasta_dict), and that the gtf transcripts
        # are all in the corresponding fasta
        self.clustered_gtf_path = clustered_gtf_path
        self.fasta_dict = fasta_dict
        self.__validate(validate_strict)

    # properties getters/setters ----------------------------------------------

    @property
    def clustered_gtf_path(self) -> str:
        """path to a gtf which has been clustered by pyranges. intended to 
        be the clustered gtf output by create_windows"""
        return self._clustered_gtf_path

    @clustered_gtf_path.setter
    def clustered_gtf_path(self, new_path: str) -> None:
        logger.debug('trying to set new clustered_gtf_path: {new_path}')
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
        """the clustered_gtf_path read in as a pr.PyRanges obj"""
        return self._clustered_gtf

    @clustered_gtf.setter
    def clustered_gtf(self, new_gtf_pr: pr.PyRanges) -> None:
        
        # validate object type
        if not isinstance(new_gtf_pr, pr.PyRanges):
            raise ValueError('clustered_gtf must be a '
                             'pyranges.PyRanges object')

        # check that unique_id is a column, if it is, ensure that it is 
        # actually unique. If it is not, add a unique_id column
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
        """the unique list of cluster levels in the gtf Cluster column"""
        return self._cluster_list

    @cluster_list.setter
    def cluster_list(self, new_list: list) -> None:
        self._cluster_list = new_list

    @property
    def fasta_dict(self) -> dict:
        """A dictionary where the keys are factor levels of the Source column 
        in the clustered_regions gtf and the value is a path to a fasta which 
        stores isoform sequences"""
        return self._fasta_dict

    @fasta_dict.setter
    def fasta_dict(self, new_fasta_dict: dict) -> None:
        """Read in the fasta paths as pysam.FastaFile objects

        Args:
            fasta_dict (dict): A dictionary where the keys are factor levels of 
            the Source column in the clustered_regions gtf and the value is 
            a path to a fasta which stores isoform sequences

        Raises:
            ValueError: raised if the fasta_dict is not a dict obj
            FileNotFoundError: raised if the fasta, __or__ the fasta.fai 
            index file does not exist
        """

        logger.debug("new dict: %s", new_fasta_dict)

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
        """read in the gtf path to a PyRanges obj

        Args:
            gtf_path (str): path to a gtf file

        Returns:
            pr.PyRanges: the gtf file read into memory as a PyRanges obj
        """
        return pr.read_gtf(gtf_path)

    def __open_fastas(self, fasta_dict) -> dict:
        """Use pysam to open the fasta file

        Args:
            fasta_dict (dict): A dictionary where the keys are factor levels of 
            the Source column in the clustered_regions gtf and the value is 
            a path to a fasta which stores isoform sequences 

        Returns:
            dict: _description_
        """
        return {k: pysam.FastaFile(v) for k, v in fasta_dict.items()}  # pylint:disable=E1101 # noqa: E501

    def __validate(self, strict: bool = True) -> None:
        """Perform validation checks on the gtf and fasta to ensure that 
        the Source and isoform/transcript names are as expected

        Args:
        strict (bool, optional): Set to False to warn, rather than 
        error, if the fasta_dict keys do not fully match the gtf Source 
        factor levels, and the transcripts in the gtf do not fully match 
        the sequence names in the fasta. Defaults to True.

        Raises:
            ValueError: Raised when an expectation on the gtf/fasta 
            correspondence is not met
        """

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

    def get_isoform_seq(self, 
                        source: str = None, 
                        isoform: str = None, 
                        unique_id: str = None, **kwargs) -> str:
        """Given either the `source` and `isoform` or the `unique_id`, 
        get the isoform sequence

        Args:
            source (str, optional): A value in the Source column of the gtf.
            if this is provided, then isoform must also be provided
            isoform (str, optional): A name of the isoform/transcript. If this 
            is provided, source must also be provided
            unique_id (str, optional): A unique_id for a given isoform. 
            if this is provided, then neither source or isoform need be.
            **kwargs (named arguments, optional): additional arguments, 
            other than reference, to pysam.FastaFile.fetch(), eg start and end

        Raises:
            IOError: Raised if the input combination is illegal, eg only 
            source or only isoform.
            KeyError: raised if source is passed and it is not a recognized 
            source
            ValueError: raised if isoform is passed and it is not a recognized 
            isoform/transcript_id

        Returns:
            str: The requested isoform sequence
        """
        
        if unique_id:
            iso_window = self.get_isoform_coord(unique_id=unique_id)
            # todo address the flake issues without ignoring
            source = iso_window.source #pylint:disable=E1101 # noqa:E261,E262
            isoform = iso_window.name #pylint:disable=E1101 # noqa:E261,E262
        
        if not source or not isoform:
            raise IOError('Either pass source and isoform, or unique_id')

        if source not in self.fasta_dict.keys():
            raise KeyError(f'the source: {source} is not in the '
                           f'fasta_dict keys: {self.fasta_dict.keys()}')
        if isoform not in set(self.clustered_gtf.transcript_id):
            raise ValueError(f'the isoform {isoform} is not in the '
                             f'transcript_id columns of clustered_gtf')

        return self.fasta_dict[source].fetch(reference=isoform, **kwargs)

    # TODO remove transcript_id field hard coding to external config
    def get_isoform_coord(self,
                          source: str = None,
                          isoform: str = None,
                          unique_id: str = None) -> Window:
        """Given either the `source` and `isoform` or the `unique_id`, 
        create a Window object which describes the isoform's coordinates

        Args:
            source (str, optional): A value in the Source column of the gtf.
            if this is provided, then isoform must also be provided
            isoform (str, optional): A name of the isoform/transcript. If this 
            is provided, source must also be provided
            unique_id (str, optional): A unique_id for a given isoform. 
            if this is provided, then neither source or isoform need be.
            **kwargs (named arguments, optional): additional arguments, 

        Raises:
            IOError: Raised if the input combination is illegal, eg only 
            source or only isoform.
            KeyError: raised if source is passed and it is not a recognized 
            source
            ValueError: raised if isoform is passed and it is not a recognized 
            isoform/transcript_id

        Returns:
            Window: _description_
        """
        # instantiate a DUMMY variable -- this will evaluate to 
        # false if the right combination of input is not found and a 
        # error will be raised
        isoform_gtf = []
        
        # The user must submit either source + isoform or unique_id. The 
        # following logic deals with that variable input. Too bad python 
        # doesn't have function overloading
        if source and isoform:
            isoform_gtf = self.clustered_gtf[
                (self.clustered_gtf.Source == source) &
                (self.clustered_gtf.transcript_id == isoform)]
        elif unique_id:
            isoform_gtf = self.clustered_gtf[
                self.clustered_gtf.unique_id == unique_id]
        else:
            raise IOError('Either source and isoform, or unique_id '
                          'must be passed in order to extract a isoform')
        # validate
        if len(isoform_gtf) == 0:
            raise ValueError('Isoform %s in Source %s not found'  #pylint:disable=C0209 # noqa:E261,E501,E262
                             % (isoform, source))
        elif len(isoform_gtf) > 1:
            raise KeyError('Isoform %s in source %s is not unique'  #pylint:disable=C0209 # noqa:E261,E501,E262
                           % (isoform, source))
        
        # convert the subsetted gtf to a dict for easy extraction
        isoform_gtf_dict = isoform_gtf.df.to_dict('records')[0]

        # create a Window object which describes the isoform
        w = Window(
            str(isoform_gtf_dict['Chromosome']),
            int(isoform_gtf_dict['Start']),
            int(isoform_gtf_dict['End']),
            str(isoform_gtf_dict['Strand']),
            name=str(isoform_gtf_dict['transcript_id']),
            source=str(isoform_gtf_dict['Source']))

        return w

    def get_cluster(self, cluster: str) -> pr.PyRanges:
        """Given a cluster id, return the cluster_gtf subsetted down 
        to just that cluster

        Args:
            cluster (str): a valid factor level in the Cluster column of the 
            gtf object

        Raises:
            ValueError: raised if the cluster id is not in the gtf Cluster 
            column

        Returns:
            pr.PyRanges: A clustered_gtf PyRanges obj for a single cluster
        """

        if cluster not in self.cluster_list:
            raise ValueError(f'{cluster} not in cluster_list')

        # NOTE: the Cluster property exists because of pyranges, though this 
        # of course depends on the pyranges API
        return self.clustered_gtf[self.clustered_gtf.Cluster == cluster]
    
    def get_cluster_coord(self, cluster: int, stranded: bool = True) -> Window:
        """Given a cluster id, create a Window object which describes 
        the cluster coordinates and number of isoforms

        Args:
            cluster (str): a valid factor level in the Cluster column of the 
            gtf object
            stranded (bool, optional): Set to False if the gtf was clustered 
            without reference to strand. Defaults to True.

        Returns:
            Window: A window object which describes a given cluster. __NOTE__ 
                    that the `score` attribute is the number of 
                    isoforms in the cluster
        """
        # extract a subset of the clustered_gtf for just 1 cluster
        cluster_gtf = self.get_cluster(cluster)
        
        # create a Window obj which describes the cluster's region
        # NOTE: the score attribute is the number of sources (individuals)
        # in the cluster. Also, this depends on the Source column being 
        # populated appropriately
        w = Window(
            list(cluster_gtf.Chromosome)[0],
            min(cluster_gtf.Start),
            max(cluster_gtf.End),
            list(cluster_gtf.Strand)[0] if stranded else "*",
            name=str(cluster),
            score=len(cluster_gtf.Source.unique()))
        
        return w
