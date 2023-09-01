# pylint:disable=C0209
import logging
import re

logger = logging.getLogger(__name__)

__all__ = ['Window']


class Window:
    """A class with properties chr (chromosome) start and end which describe 
    a 0 based half open interval, ie [0,100) describes a window of length 100 
    with the first accessed item at index 0 and the last accessed item at 
    index 99. This mimics the BED format standard."""

    _STRAND_VALUES = ['+', '-', '*']

    def __init__(self, chr: str, window_start: int,
                 window_end: int, strand: str,
                 name: str = "", score: int = 1000, **kwargs):
        self._chr = chr
        self._start = window_start
        self._end = window_end
        self._strand = strand
        self._name = name
        self._score = score
        
        for k, v in kwargs.items():
            setattr(self, k, v)

    def __len__(self):
        """length of the region"""
        return self._end - self._start
    
    # TODO print additional keys which may be passed after this bed6 formatted 
    # string
    def __str__(self):
        """print a description of the region to std out"""
        return "\t".join([self._chr, str(self._start),
                         str(self._end), str(self.name), str(self.score),
                          self._strand])

    @property
    def chr(self):
        """The chromosome. This either includes a prefix 'chr' or not"""
        return self._chr

    @chr.setter
    def chr(self, new_chr: str):
        if not isinstance(new_chr, str):
            logger.debug(new_chr)
            raise ValueError('chr must be a string')
        self._chr = new_chr

    @property
    def start(self):
        """This is the coorindate of where our bin/window begins. Integer. It 
        is a 1-based index location on the chromosome."""
        return self._start

    @start.setter
    def start(self, new_start: int):
        if not isinstance(new_start, int):
            logger.debug(new_start)
            raise ValueError('start must be an integer')
        self._start = new_start

    @property
    def end(self):
        """This is the coorindate of where our bin/window begins. Integer. 
        It is a 1-based index location on the chromosome."""
        return self._end

    @end.setter
    def end(self, new_end: int):
        if not isinstance(new_end, int):
            logger.debug(new_end)
            raise ValueError('end must be an integer')
        self._end = new_end

    @property
    def strand(self):
        """This is the coorindate of where our bin/window begins. Integer. 
        It is a 1-based index location on the chromosome."""
        return self._strand

    @strand.setter
    def strand(self, new_strand: str):
        if new_strand not in self._STRAND_VALUES:
            logger.debug(new_strand)
            raise ValueError('strand value: %s is not one of the recognized '
                             'strand values: %s'
                             % (new_strand, ','.join(self._STRAND_VALUES)))
        self._strand = new_strand

    @property
    def name(self):
        """Store the name of a given window. Default constructor sets this to 
        an empty string"""
        return self._name

    @name.setter
    def name(self, new_name: int):
        new_name = str(new_name)
        self._name = new_name

    @property
    def score(self):
        """Store a score for a given window. Default constructor sets this 
        to 1000"""
        return self._score

    @score.setter
    def score(self, new_score: int):

        logger.debug('new score: %s', new_score)

        if not isinstance(new_score, int):
            raise ValueError('Score must be an integer')
        elif new_score < 0 or new_score > 1000:
            Warning('Bed format requires score in interval [0,1000]')

        self._score = new_score
    
    def to_dict(self) -> dict:
        """transform a window's attributes into a dictionary

        Returns:
            dict: The attributes of a given window reformatted as a dictionary
        """
        d = vars(self)
        # remove leading underscore from keys
        return {re.sub(r"^_", '', k): v for k, v in d.items()}
