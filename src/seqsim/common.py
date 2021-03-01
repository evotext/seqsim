"""
Common functions.
"""

# Import Python standard libraries
import hashlib
import random
from typing import Hashable, Union, List, Optional, Sequence

# Import 3rd party libraries
import numpy as np


def set_seeds(seed: Union[str, float, int]) -> None:
    """
    Set seeds globally from the one provided by the user.

    The function takes care of reproducibility and allows to use strings and
    floats as seed for `numpy` as well.

    :param seed: The seed for Python and numpy random number generators.
    """

    # Set seed for Python RNG
    random.seed(seed)

    # Allows using strings as numpy seeds, which only takes uint32 or arrays of uint32
    if isinstance(seed, (str, float)):
        np_seed = np.frombuffer(
            hashlib.sha256(str(seed).encode("utf-8")).digest(), dtype=np.uint32
        )
    else:
        np_seed = seed

    # Set the np set
    np.random.seed(np_seed)


# TODO: properly rewrite
def sequence_find(hay: Sequence, needle: Sequence) -> Optional[int]:
    """
    Return the index for starting index of a sub-sequence within a sequence.

    The function is intended to work similarly to the built-in `.find()` method for
    Python strings, but accepting all types of sequences (including different types
    for `hay` and `needle`).

    :param hay: The sequence to be searched within.
    :param needle: The sub-sequence to be located in the sequence.
    :return: The starting index of the sub-sequence in the sequence, or `None` if not
             found.
    """
    # Cache `needle` length and have it as a tuple already
    len_needle = len(needle)
    t_needle = tuple(needle)

    # Iterate over all sub-lists (or sub-tuples) of the correct length and check
    # for matches
    for i in range(len(hay) - len_needle + 1):
        if tuple(hay[i : i + len_needle]) == t_needle:
            return i

    return None


# TODO: replace with the ngram collector module
def collect_subseqs(sequence: Sequence, sort: bool = True) -> List[Sequence]:
    """
    Collects all possible sub-sequences in a given sequence.

    When sorting is requested, sub-sequences will first be sorted by their length and,
    later, by comparing one with the other. Mixing types, like strings and integers, can
    lead to unexpected results and is not suggested if the type cannot be guaranteed.

    Note that this function performs simple comprehensions, neither using padding
    symbols nor the more complex methods n-gram collection methods ultimately based on
    `ngram_iter()`.

    Examples
    --------
    collect_subseqs('abcde')
    ['a', 'b', 'c', 'd', 'e', 'ab', 'bc', 'cd', 'de', 'abc', 'bcd', 'cde', 'abcd', 'bcde', 'abcde']

    :param sequence: The sequence that shall be converted into it's ngram-representation.
    :param sort: Whether to sort the list of ngrams by length and by identity
        (default: True).
    :return: A list of all ngrams of the input sequence.
    """

    # Cache the length of the sequence
    length = len(sequence)

    # Set the starting index
    idx = 0

    # define the output list
    ret = []

    # start the while loop
    while idx != length and idx < length:
        # copy the sequence
        new_sequence = sequence[idx:length]

        # append the sequence to the output list
        ret += [new_sequence]

        # loop over the new sequence
        for j in range(1, len(new_sequence)):
            ret += [new_sequence[:j]]
            ret += [new_sequence[j:]]

        # increment idx and decrement length
        idx += 1
        length -= 1

    if sort:
        # We try to sort normally; if there is a TypeError, such as when the list has mixed
        # ints and strings, we sort by the string representation of all elements
        # TODO: do it in a better way
        try:
            ret = sorted(ret, key=lambda e: (len(e), e))
        except TypeError:
            ret = sorted(ret, key=lambda e: (len(str(e)), str(e)))

    return ret


def _nwise(I, n):
    """
    Iterate through I n at a time, e.g.
        _nwise("abcd", 2) -> ab, bc, cd
    """
    for i in range(len(I) + 1 - n):
        yield I[i : i + n]


def _indices(L: Sequence[Hashable], element: Hashable) -> list:
    """
    Find all _indices in `L` matching `element`, e.g.
        _indices("abcab", "a") -> 0, 3
    """

    return [idx for idx, value in enumerate(L) if value == element]
