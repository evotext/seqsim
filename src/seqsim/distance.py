"""
Module implementing various methods for distance in sequences.

Most of the methods are commonly used in string comparison, such as Levenshtein
distance, but in this module we need to make sure we can operate on arbitrary
iterable data structures.
"""

# TODO: redo distance/similarity

# Import Python standard libraries
from typing import Sequence, List, Tuple

# Import 3rd party libraries
import numpy as np

# Import local modules
from .common import sequence_find, collect_subseqs


# TODO: allow some normalization?
# TODO: custom costs?
def edit_distance(seq1: Sequence, seq2: Sequence) -> float:
    """
    Computes a common Levenshtein edit distance between two sequences.

    :param seq1: The first sequence of elements to be compared.
    :param seq2: The second sequence of elements to be compared.
    :return: The Levenshtein distance between the two sequences.
    """

    # For all `x` and `y`, matrix[i, j] will hold the Levenshtein distance between the
    # first `x` characters of `seq1` and the first `y` characters of `seq2`; the
    # starting matrix is a zeroed one
    size_x = len(seq1) + 1
    size_y = len(seq2) + 1
    matrix = np.zeros((size_x, size_y))

    # Source prefixes can be transformed into empty sequences by dropping all chars
    for x in range(size_x):
        matrix[x, 0] = x

    # Target prefixes can be reached from empty sequence prefix by inserting every char
    for y in range(size_y):
        matrix[0, y] = y

    # Main loop, with the implied substitution cost
    for x in range(1, size_x):
        for y in range(1, size_y):
            if seq1[x - 1] == seq2[y - 1]:
                matrix[x, y] = min(
                    matrix[x - 1, y] + 1,  # deletion
                    matrix[x - 1, y - 1],  # insertion
                    matrix[x, y - 1] + 1,  # substitution
                )
            else:
                matrix[x, y] = min(
                    matrix[x - 1, y] + 1,  # deletion
                    matrix[x - 1, y - 1] + 1,  # insertion
                    matrix[x, y - 1] + 1,  # substitution
                )

    return matrix[size_x - 1, size_y - 1]


def jaccard_distance(seq1: Sequence, seq2: Sequence) -> float:
    """
    Computes a Jaccard distance between two sequences.

    :param seq1: The first sequence of elements to be compared.
    :param seq2: The second sequence of elements to be compared.
    :return: The Jaccard distance between the two sequences.
    """

    intersection = len(set(seq1).intersection(seq2))
    union = len(seq1) + len(seq2) - intersection
    return 1.0 - (float(intersection) / union)


# TODO: rename from jaccard, as it is not really jaccard...
def subseq_jaccard_distance(seq1: Sequence, seq2: Sequence) -> float:
    """
    Computes a Jaccard distance between two sequences using sub-sequence occurrence.

    :param seq1: The first sequence of elements to be compared.
    :param seq2: The second sequence of elements to be compared.
    :return: The Subseq-Jaccard distance between the two sequences.
    """

    subseqs1 = collect_subseqs(seq1)
    subseqs2 = collect_subseqs(seq2)

    # From the longest subseq, which is the length of the longest sequence,
    # collect all subsequences of that given length in both sets, compute the
    # Jaccard index, correct it by length of the subsequence (so that longer ones
    # will score higher) and update the internal results before returning.
    # TODO: collect beforehand in a list, so we don´t repeat the comprehension
    # TODO: convert to tuple beforehand as well?
    jaccard = []
    max_length = max([len(seq1), len(seq2)])
    for length in range(max_length, 0, -1):
        l_subseq1 = [tuple(ss) for ss in subseqs1 if len(ss) == length]
        l_subseq2 = [tuple(ss) for ss in subseqs2 if len(ss) == length]

        intersection = len(set(l_subseq1).intersection(l_subseq2))
        union = len(l_subseq1) + len(l_subseq2) - intersection
        jaccard.append((float(intersection) / union) * length)

    # Compute the denominator, as the highest possible value
    den = (max_length * (max_length + 1)) / 2.0

    return (1.0 - (sum(jaccard) / den)) ** max_length


def _mmcwpa(
    f_x: List[Sequence], f_y: List[Sequence], ssnc: float
) -> Tuple[List[Sequence], List[Sequence], float]:
    """
    Internal function for MMCWPA implementation.

    In this implementation of the Modified Moving Contracting Window Pattern Algorithm
    (MMCWPA) to calculate sequence similarity, we return a list of non-overlapping,
    non-contiguous fields Fx, a list of non-overlapping, non-contiguous fields Fy, and
    the SSNC value (the Sum of the Square of the Number of the same characters). This
    function separates the core method of the implementation and makes recursive calls
    easier.

    :param f_x: A list of sub-sequences, related to the first sequence.
    :param f_y: A list of sub-sequences, related to the second sequence.
    :param ssnc: The previous SSNC value.
    :return: A tuple whose first element is a list of remaining sub-sequences from the
             first sequence, the second element is a list of remaining sub-sequences
             from the second sequence, and the third element is the updated SSNC.
    """

    # the boolean value indicating if a total or partial
    # match was found between subfields Fx and Fy; when
    # a match is found, the variable is used to cascade
    # out of the loops of the function
    match = False

    # the variables where to store the new collections of
    # subfields, if any match is found; if these values
    # are not changed and the empty lists are returned,
    # stringcomp() will break the loop of comparison,
    # calculate the similarity ratio and return its value
    new_f_x, new_f_y = [], []

    # search patterns in all subfields of Fx; the index of
    # the subfield in the list is used for upgrading the
    # list, if a pattern is a found
    for idx_x, sf_x in enumerate(f_x):
        # 'length' stores the length of the sliding window,
        # from full length to a single character
        for length in range(len(sf_x), 0, -1):
            # 'i' stores the starting index of the sliding
            # window in Fx
            for i in range(len(sf_x) - length + 1):
                # extract the pattern for matching
                pattern = sf_x[i : i + length]

                # look for the pattern in Fy
                for idx_y, sf_y in enumerate(f_y):
                    # 'j' stores the starting index in Fy; the
                    # Python find() function returns -1 if there
                    # is no match
                    j = sequence_find(sf_y, pattern)
                    if j is not None:
                        # the pattern was found; set 'new_fx' and
                        # 'new_fy' to version of 'fx' and 'fy' with
                        # the patterns removed, update the SSNC and
                        # set 'match' as True, in order to cascade
                        # out of the loops
                        tmp_x = [sf_x[:i], sf_x[i + length :]]
                        tmp_y = [sf_y[:j], sf_y[j + length :]]
                        new_f_x = f_x[:idx_x] + tmp_x + f_x[idx_x + 1 :]
                        new_f_y = f_y[:idx_y] + tmp_y + f_y[idx_y + 1 :]

                        ssnc += (2 * length) ** 2

                        match = True

                        break

                    # if the current match was found, end search
                    if match:
                        break

                # if a match was found, end the sliding window
                if match:
                    break

            # if a match was found, end Fx subfield enumeration
            if match:
                break

        # remove any empty subfields due to pattern removal
        new_f_x = [sf for sf in new_f_x if sf]
        new_f_y = [sf for sf in new_f_y if sf]

        return new_f_x, new_f_y, ssnc


def mmcwpa_distance(seq_a: Sequence, seq_b: Sequence) -> float:
    """
    Computes an MMCWPA distance between two sequences.

    MMCWPA is the Modifier Moving Contracting Window Pattern Algorithm, modified by
    Tiago Tresoldi from a method published by Yang et al. (2001). In order to simplify
    the logic, the function uses the auxiliary internal function `_mmcwpa()`.

    Reference for the original method:
    Q. X. Yang, S. S. Yuan, L. Chun, L. Zhao, S. Peng. "Faster Algorithm of String
    Comparison", eprint arXiv:cs/0112022, December 2001.

    :param seq_a: The first sequence of elements to be compared.
    :param seq_b: The second sequence of elements to be compared.
    :return: The MMCWPA distance between the two sequences.

    """

    len_x, len_y = len(seq_a), len(seq_b)

    f_x, f_y = [seq_a], [seq_b]

    ssnc = 0.0
    while True:
        f_x, f_y, ssnc = _mmcwpa(f_x, f_y, ssnc)

        if len(f_x) == 0 or len(f_y) == 0:
            break

    return 1.0 - ((ssnc / ((len_x + len_y) ** 2.0)) ** 0.5)
