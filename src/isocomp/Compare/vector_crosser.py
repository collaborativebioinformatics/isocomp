import operator
import logging
from itertools import product

import numpy as np

logger = logging.getLogger(__name__)

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


# TODO check the two functions below for equivalent functionality. If 
# they are the same, then replace the old implementation with one of the 
# new ones -- preferrably the stdlib unless there is a compelling reason not 
# to. Consider returning a list of tuples instead of the dict.

# this needs to be checked...if it achieves the same thing, then replace
# the old implementation. Returning the combinations (list of tuples),
# rather than the dict, is also better
def __vector_crosser_stdlib(v1: list, v2: list, equals: bool = False) -> dict:
    """
    Given two lists with any length and any element type, generate a
    dictionary with keys 'V1' and 'V2', each of which stores a list.
    Indices of the list correspond to one another which describe all
    unique combinations of the elements of v1 and v2.
    Set equals to TRUE to return corresponding elements with equal values,
    e.g., 1 == 1. This is based on R code here:
    https://github.com/mhesselbarth/suppoRt/blob/HEAD/R/expand_grid_unique.R

    Args:
        v1 (list): a list of items
        v2 (list): a list of items
        equals (bool, optional): whether to return paired elements where
        the values of v1 and v2 are the same, e.g., '1' '1' would be in the 
        same index in V1 and V2 if this is set to True. Defaults to False.

    Returns:
        dict: a dictionary with keys 'V1' and 'V2', each of which stores a
        list. Indices of the list correspond to one another which describe
        all unique combinations of the elements of v1 and v2
    """
    unique_v1 = list(set(v1))
    unique_v2 = list(set(v2))

    if equals:
        combinations = list(product(unique_v1, unique_v2))
    else:
        combinations = [(a, b) for a in unique_v1 for b in unique_v2 if a != b]

    d = {
        "V1": [x[0] for x in combinations],
        "V2": [x[1] for x in combinations]
    }

    return d


def __vector_crosser_numpy(v1: list, v2: list, equals: bool = False) -> dict:
    """
    Given two lists with any length and any element type, generate a
    dictionary with keys 'V1' and 'V2', each of which stores a list.
    Indices of the list correspond to one another which describe all
    unique combinations of the elements of v1 and v2.
    Set equals to TRUE to return corresponding elements with equal values,
    e.g., 1 == 1. This is based on R code here:
    https://github.com/mhesselbarth/suppoRt/blob/HEAD/R/expand_grid_unique.R

    Args:
        v1 (list): a list of items
        v2 (list): a list of items
        equals (bool, optional): whether to return paired elements where
        the values of v1 and v2 are the same, e.g., '1' '1' would be in the same
        index in V1 and V2 if this is set to True. Defaults to False.

    Returns:
        dict: a dictionary with keys 'V1' and 'V2', each of which stores a
        list. Indices of the list correspond to one another which describe
        all unique combinations of the elements of v1 and v2
    """
    unique_v1 = np.array(list(set(v1)))
    unique_v2 = np.array(list(set(v2)))

    v1_grid, v2_grid = np.meshgrid(unique_v1, unique_v2, indexing='ij')
    combinations = np.column_stack((v1_grid.ravel(), v2_grid.ravel()))

    if not equals:
        mask = combinations[:, 0] != combinations[:, 1]
        combinations = combinations[mask]

    d = {
        "V1": list(combinations[:, 0]),
        "V2": list(combinations[:, 1])
    }

    return d