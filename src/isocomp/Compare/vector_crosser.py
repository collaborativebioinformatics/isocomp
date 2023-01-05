import operator
import logging

from .IsoformLibrary import IsoformLibrary

logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = ['vector_crosser']


# TODO this isn't very efficient b/c of the list operations. There is likely a
# better implementation maybe in numpy, or a more pythonic way of doing this
# same thing in a lot less code
def vector_crosser(v1: list, v2: list, equals: bool = False) -> dict:
    """given two lists with any length and any element type, generate a  
    a dictionary with keys 'V1' and 'V2', each of which stores a list. 
    Indicies of the list correspond to one another which describe all 
    unique combinations of the elements of v1 and v2. 
    Set equals to TRUE to return corresponding elements with equal values, 
    eg 1 == 1. This is based on R code here:
    https://github.com/mhesselbarth/suppoRt/blob/HEAD/R/expand_grid_unique.R

    Args:
        v1 (list): a list of items
        v2 (list): a list of items
        equals (bool, optional): whether to return paired elements where 
        the values of v1 and v2 are the same, eg '1' '1' would be in the same 
        index in V1 and V2 if this is set to True. Defaults to False.

    Returns:
        dict: a dictionary with keys 'V1' and 'V2', each of which stores a 
        list. Indicies of the list correspond to one another which describe 
        all unique combinations of the elements of v1 and v2
    """
    d = {}

    unique_v1 = list(set(v1))
    unique_v2 = list(set(v2))

    def inner(i: int) -> None:
        """This is intended to be used in the for loop below. The variable 
        z stores the set diff between unique_v2 and, depending on the value 
        of i and the variable equals, some range of unique_v1. For example, 
        in the for loop below, we iterate over the length of unique_v1. If the 
        length is three, __and__ equals is set to False, then the first 
        iteration takes the set diff of unique_v2 and unique_v1[0:1] which 
        is the first element of unique_v1. If equals is set to True, then the 
        first iteration is the set diff of unique_v2 and unique_v1[0:0] which 
        returns the entirety of unique_v2. this continues in the for loop below, 
        iteratively taking more of unique_v1

        Args:
            i (int): This is used to extract a range of unique_v1 in the 
            set difference operation, and to extract a a given value from 
            unique_v1 and append it (repeated for length(z)) to V1 while 
            z (the set diff result) is appended to V2
        """
        z = list(set(unique_v2) - set(unique_v1[0:i + operator.not_(equals)]))
        if z:
            d.setdefault('V1', []).extend([unique_v1[i]]*len(z))
            d.setdefault('V2', []).extend(z)

    # see the docstring for inner() above
    for i in range(len(unique_v1)):
        inner(i)

    return d
