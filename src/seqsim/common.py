"""
Common functions.
"""

# Import Python standard libraries
import hashlib
import random
from typing import Callable, Hashable, Union, List, Optional, Sequence, Tuple
import string

# Import 3rd party libraries
import numpy as np

# TODO: support groups with more than 100 elements (length of string.printable)
def equivalent_string(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]
) -> Tuple[str, str]:
    """
    Returns a string equivalent to a sequence, for comparison.

    As some methods offered by third-party libraries only operate on
    strings, while `seqsim` is designed to offer all methods of
    comparison for generic sequences of hashable elements, in
    some cases it is necessary to convert a sequence to an equivalent
    string. Using a normal `str` conversion is not possible or
    satisfactory for a number of reasons, including elements not
    having a string representation, and individual string
    representations of different lengths and potentially overlapping
    (consider cases like `[1, 12, 123, 23]`).

    This function accepts a pair of sequences and returns an equivalent
    textual representation, that is, a pair of strings where the
    order is preserved and each token is mapped to a single, unique
    character. While the information in the strings is meaningless,
    they are built to facilitate inspection and debugging as much
    as possible, trying to use only ASCII printable characters or
    Unicode characters that are expected to be supported for
    visualization in the majority of systems.

    If two strings are passed, the same strings will be returned. Note
    that in case of mixed types (such as a string and a list of
    strings), strings will be considered sequences of characters
    and will be modified upon return, as the tokens of the
    second sequence could be of length over one character (e.g.,
    `"abc"` and `["a", "bc"]`).

    :param seq_x: The first sequence to be mapped to an equivalent
        string.
    :param seq_y: The second sequence to be mapped to an equivalent
        string.
    :return: A tuple of two strings equivalent, for matters of
        comparison and distance computation, to the provided
        sequences.
    """

    # Don't need to apply to strings
    if isinstance(seq_x, str) and isinstance(seq_y, str):
        return seq_x, seq_y

    # Map the sequences to lists and get the set of symbols
    # that are used (the list is sorted for reproducibility)
    seq_x = [element for element in seq_x]
    seq_y = [element for element in seq_y]
    elements = sorted(set(seq_x + seq_y), key=lambda e: str(e))

    # Use ASCII printable elements if possible
    if len(elements) <= len(string.printable):
        mapper = {source: target for source, target in zip(elements, string.printable)}
    else:
        raise ValueError("Too many values")

    # Map the new sequences with string elements and return
    new_seq_x = [mapper.get(element) for element in seq_x]
    new_seq_y = [mapper.get(element) for element in seq_y]

    return "".join(new_seq_x), "".join(new_seq_y)


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


# TODO: continue `Callable` typing
def _wagner_fischer(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    costs_fn: Callable,
    d: Optional[List[List[float]]] = None,
) -> float:
    """
    Computes arbitrary edit distances using the Wagner-Fischer algorithm.

    See: https://en.wikipedia.org/wiki/Wagner%E2%80%93Fischer_algorithm

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param costs_fn: A cost function, specific to a particular edit distance.
    :param d: An optional "starting matrix".
    :return: The cost distance.
    """

    # Cache lengths
    m = len(seq_x)
    n = len(seq_y)

    # If we haven't been provided a custom initial matrix specific to
    # our particular distance measure, create a generic on here.  This
    # fills out the first row and first column of the cost matrix,
    # corresponding to inserting or deleting all characters.
    if not d:
        d = [[0 for i in range(0, n + 1)] for j in range(0, m + 1)]
        for i in range(1, m + 1):
            d[i][0] = i
        for j in range(1, n + 1):
            d[0][j] = j

    # Flood fill the rest of the matrix with values computed using our
    # cost function
    for j in range(1, n + 1):
        for i in range(1, m + 1):
            d[i][j] = min(costs_fn(seq_x, seq_y, d, i, j))

    # Return lower-right element of matrix, which is the minimum cost
    # for transforming s into t
    return d[m][n]
