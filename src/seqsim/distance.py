"""
Module implementing various methods for distance in sequences.

Most of the methods are commonly used in string comparison, such as Levenshtein
distance, but in this module we need to make sure we can operate on arbitrary
iterable data structures.

We follow the mathematical definitions for distinguishing between "similarity"
and "distance", as the latter must have the following properties:

  * positivity: d(x,y) >= 0
  * symmetry: d(x,y) = d(y,x)
  * identity-discerning: d(x,y) = 0 => x = y
  * triangle inequality: d(x,z) <= d(x,y) + d(y,z)
"""

# TODO: allow some normalization for all methods?

# Import Python standard libraries
from typing import Hashable, List, Sequence, Tuple

# Import 3rd party libraries
import numpy as np

# Import local modules
from .common import sequence_find, collect_subseqs

from . import similarity

# TODO: review code, especially (N,N) (shoud probably delegate it to similarity)
def fast_birnbaum(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
    if len(seq_y) > len(seq_x):
        seq_x, seq_y = seq_y, seq_x

    distance = 1.0 - (
        similarity.fast_birnbaum_simil(seq_x, seq_y)
        / similarity.fast_birnbaum_simil(seq_y, seq_y)
    )

    # Distance can be less than zero if M contains a perfect N plus doublets
    # matching part of N as well.
    return max(0.0, distance)


def levenshtein(seq_a, seq_b):
    """
    Compute the Levenshtein distance between two manuscripts, using the
    Wagner-Fischer algorithm.
    See: https://en.wikipedia.org/wiki/Levenshtein_distance
    """
    return _wagner_fischer(seq_a, seq_b, _levenshtein_costs)


def _levenshtein_costs(seq_a, seq_b, d, i, j):
    """
    Compute the candidate costs for entry (i,j) of the Levenshtein distance
    matrix, each cost corresponding to one of the available edit operations.
    """
    substitution_cost = 0 if seq_a[i - 1] == seq_b[j - 1] else 1
    costs = (
        d[i - 1][j] + 1,  # del
        d[i][j - 1] + 1,  # ins
        d[i - 1][j - 1] + substitution_cost,
    )
    return costs


def norm_levenshtein(seq_a, seq_b):
    """
    Compute a normalised form of the Levenshtein distance between two
    manuscripts, guaranteeing a result between 0 and 1 regardless of
    manuscript length.
    """
    return levenshtein(seq_a, seq_b) / max([len(seq_a), len(seq_b)])


def levdamerau(seq_a, seq_b):
    """
    Compute the Damerau-Levenshtein distance between two manuscripts,
    using the Wagner-Fischer algorithm.
    See: https://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance
    """
    return _wagner_fischer(seq_a, seq_b, _levdamerau_costs)


def _levdamerau_costs(seq_a, seq_b, d, i, j):
    """
    Compute the candidate costs for entry (i,j) of the Damerau-Levenshtein
    distance matrix, each cost corresponding to one of the available edit
    operations.
    """
    # Start out as per Levenshtein
    substitution_cost = 0 if seq_a[i - 1] == seq_b[j - 1] else 1
    costs = [
        d[i - 1][j] + 1,
        d[i][j - 1] + 1,
        d[i - 1][j - 1] + substitution_cost,
    ]
    # Add the transposition cost
    if (
        i > 1
        and j > 1
        and seq_a[i - 1] == seq_b[j - 2]
        and seq_a[i - 2] == seq_b[j - 1]
    ):
        costs.append(d[i - 2][j - 2] + 1)
    return costs


def fragile_ends(seq_a, seq_b):
    """
    Compute the "fragile ends" distance between two manuscripts - like
    Levenshtein but with deletions in the initial or final 10% of the
    positions being discounted.
    This distance measure is not used directly in the paper and was a
    proof-of-concept while working toward the "Stemmatological distance".
    """
    d = _fragile_ends_initial_matrix(seq_a, seq_b)
    return _wagner_fischer(seq_a, seq_b, _fragile_ends_costs, d)


def _fragile_ends_initial_matrix(seq_a, seq_b):
    """
    Compute the starting matrix for the Wagner-Fischer algorithm when using
    "fragile ends" costs.  The generic starting matrix isn't applicable here,
    as the cost of deleting the entire starting manuscript is less than the
    manuscript length.
    """
    m = len(seq_a)
    n = len(seq_b)
    d = [[0 for i in range(0, n + 1)] for j in range(0, m + 1)]
    for i in range(1, m + 1):
        if i <= round(0.1 * m) or i >= round(0.9 * m):
            d[i][0] = d[i - 1][0] + 0.5
        else:
            d[i][0] = d[i - 1][0] + 1
    # Insertions
    for j in range(1, n + 1):
        d[0][j] = j
    return d


def _fragile_ends_costs(seq_a, seq_b, d, i, j):
    """
    Compute the candidate costs for entry (i,j) of the "fragile ends"
    distance matrix, each cost corresponding to one of the available edit
    operations.
    """
    substitution_cost = 0 if seq_a[i - 1] == seq_b[j - 1] else 1
    costs = [
        d[i][j - 1] + 1,
        d[i - 1][j - 1] + substitution_cost,
    ]
    # Discount deletion near ends
    m = len(seq_a)
    if i <= round(0.1 * m) or i >= round(0.9 * m):
        costs.append(d[i - 1][j] + 0.5)
    else:
        costs.append(d[i - 1][j] + 1)
    return costs


def bulk_delete(seq_a, seq_b, max_del_len=5):
    """
    Compute the "bulk delete" distance between two manuscripts - like
    Levenshtein but with deletions of consecutive stories up to the specified
    maximum length counting as a single operation.
    This distance measure is not used directly in the paper and was a
    proof-of-concept while working toward the "Stemmatological distance".
    """
    d = _bulk_delete_initial_matrix(seq_a, seq_b, max_del_len)
    _bulk_delete_costs = _bulk_delete_costs_factory(max_del_len)
    return _wagner_fischer(seq_a, seq_b, _bulk_delete_costs, d)


def _bulk_delete_initial_matrix(seq_a, seq_b, max_del_len):
    """
    Compute the starting matrix for the Wagner-Fischer algorithm when using the
    "bulk delete distance" costs.  The generic starting matrix isn't applicable
    here, as the cost of deleting the entire starting manuscript is less than
    the manuscript length.
    """
    m = len(seq_a)
    n = len(seq_b)
    d = [[0 for i in range(0, n + 1)] for j in range(0, m + 1)]
    # Bulk deletions
    for i in range(1, m + 1):
        bulk_deletions = int(i / max_del_len)
        d[i][0] = bulk_deletions
        remainder = i - bulk_deletions * max_del_len
        if remainder:
            d[i][0] += 1
    # Insertions
    for j in range(1, n + 1):
        d[0][j] = j
    return d


def _bulk_delete_costs_factory(max_del_len=5):
    """
    Define and return a function which will compute the candidate costs for
    entry (i,j) of the "bulk deletion" distance matrix with a particular maximum
    deletion block length.
    """

    def _bulk_delete_costs(seq_a, seq_b, d, i, j):
        substitution_cost = 0 if seq_a[i - 1] == seq_b[j - 1] else 1
        costs = [
            d[i][j - 1] + 1,  # ins
            d[i - 1][j - 1] + substitution_cost,
        ]
        for n in range(1, min(max_del_len + 1, i)):
            # Delete consecutive block of n
            costs.append(d[i - n][j] + 1)
        return costs

    return _bulk_delete_costs


def stemmatological(seq_a, seq_b, max_del_len=5, frag_start=10, frag_end=10):
    """
    Compute the "stemmatological distance" between two manuscripts.  This is
    essentially a combination of the above "fragile ends" and "bulk delete"
    methods, with the fragile ends method generalised a little to allow
    specifying the size of both fragile regions.
    """
    d = _stemmatological_initial_matrix(seq_a, seq_b, max_del_len, frag_start, frag_end)
    _stemmatology_costs = _stemmatological_costs_factory(
        max_del_len, frag_start, frag_end
    )
    return _wagner_fischer(seq_a, seq_b, _stemmatology_costs, d)


def _stemmatological_initial_matrix(
    seq_a, seq_b, max_del_len=5, frag_start=10, frag_end=10
):
    """
    Compute the starting matrix for the Wagner-Fischer algorithm when using the
    "stemmatological distance" costs.  The generic starting matrix isn't
    applicable here, as per the fragile ends and bulk delete methods above.
    """
    m = len(seq_a)
    n = len(seq_b)
    d = [[0 for i in range(0, n + 1)] for j in range(0, m + 1)]
    lower = round(m * (frag_start) / 100.0)
    upper = round(m * (100 - frag_end) / 100.0)
    for i in range(1, m + 1):
        if i <= lower or i >= upper:
            d[i][0] = d[i - min(i, max_del_len)][0] + 0.5
        else:
            d[i][0] = d[i - min(i, max_del_len)][0] + 1.0
    # Insertions
    for j in range(1, n + 1):
        d[0][j] = j
    return d


def _stemmatological_costs_factory(max_del_len=5, frag_start=10, frag_end=10):
    """
    Define and return a function which will compute the candidate costs for
    entry (i,j) of the "stemmatological distance" matrix with specific parameter
    values.
    """

    def _stemmatological_costs(seq_a, seq_b, d, i, j):
        substitution_cost = 0 if seq_a[i - 1] == seq_b[j - 1] else 1
        costs = [
            d[i][j - 1] + 1,
            d[i - 1][j - 1] + substitution_cost,
        ]
        m = len(seq_a)
        lower = round(m * (frag_start) / 100.0)
        upper = round(m * (100 - frag_end) / 100.0)
        # Delete consecutive block of n
        for n in range(1, min(max_del_len, i)):
            # Discount bulk deletion near ends
            if i <= lower or i >= upper:
                costs.append(d[i - n][j] + 0.5)
            else:
                costs.append(d[i - n][j] + 1)
        return costs

    return _stemmatological_costs


def norm_stemmatological(seq_a, seq_b):
    """
    Compute a normalised form of the "stemmatological distance"
    between two manuscripts, guaranteeing a result between 0 and 1 regardless of
    manuscript length.
    """
    return stemmatological(seq_a, seq_b) / max([len(seq_a), len(seq_b)])


def norm_stemmatological_2030(seq_a, seq_b):
    """
    Compute a normalised form of the "stemmatological distance"
    between two manuscripts, with the first 20% and final 30% of the manuscript
    being "fragile".
    """
    return stemmatological(seq_a, seq_b, 5, 20, 30) / max([len(seq_a), len(seq_b)])


def _wagner_fischer(seq_a, seq_b, costs_fn, d=None):
    """
    Implements the Wagner-Fischer algorithm to compute arbitrary edit
    distances.
    s is the "source" sequence (of length m)
    t is the "target" sequence (of length n)
    costs_fn is a cost function, specific to a particular edit distance
    d is an optional "starting matrix"
    See: https://en.wikipedia.org/wiki/Wagner%E2%80%93Fischer_algorithm
    """
    m = len(seq_a)
    n = len(seq_b)

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
            costs = costs_fn(seq_a, seq_b, d, i, j)
            d[i][j] = min(costs)

    # Return lower-right element of matrix, which is the minimum cost
    # for transforming s into t
    return d[m][n]


def jaccard(seq_a: Sequence, seq_b: Sequence) -> float:
    """
    Computes a Jaccard distance between two sequences.

    :param seq_a: The first sequence of elements to be compared.
    :param seq_b: The second sequence of elements to be compared.
    :return: The Jaccard distance between the two sequences.
    """

    intersection = len(set(seq_a).intersection(seq_b))
    union = len(seq_a) + len(seq_b) - intersection
    return 1.0 - (float(intersection) / union)


# TODO: rename from jaccard, as it is not really jaccard...
def subseq_jaccard(seq_a: Sequence, seq_b: Sequence) -> float:
    """
    Computes a Jaccard distance between two sequences using sub-sequence occurrence.

    :param seq_a: The first sequence of elements to be compared.
    :param seq_b: The second sequence of elements to be compared.
    :return: The Subseq-Jaccard distance between the two sequences.
    """

    subseqs1 = collect_subseqs(seq_a)
    subseqs2 = collect_subseqs(seq_b)

    # From the longest subseq, which is the length of the longest sequence,
    # collect all subsequences of that given length in both sets, compute the
    # Jaccard index, correct it by length of the subsequence (so that longer ones
    # will score higher) and update the internal results before returning.
    # TODO: collect beforehand in a list, so we don´t repeat the comprehension
    # TODO: convert to tuple beforehand as well?
    jaccard = []
    max_length = max([len(seq_a), len(seq_b)])
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
    f_a: List[Sequence], f_b: List[Sequence], ssnc: float
) -> Tuple[List[Sequence], List[Sequence], float]:
    """
    Internal function for MMCWPA implementation.

    In this implementation of the Modified Moving Contracting Window Pattern Algorithm
    (MMCWPA) to calculate sequence similarity, we return a list of non-overlapping,
    non-contiguous fields Fx, a list of non-overlapping, non-contiguous fields Fy, and
    the SSNC value (the Sum of the Square of the Number of the same characters). This
    function separates the core method of the implementation and makes recursive calls
    easier.

    :param f_a: A list of sub-sequences, related to the first sequence.
    :param f_b: A list of sub-sequences, related to the second sequence.
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
    new_f_a, new_f_b = [], []

    # search patterns in all subfields of Fx; the index of
    # the subfield in the list is used for upgrading the
    # list, if a pattern is a found
    for idx_a, sf_a in enumerate(f_a):
        # 'length' stores the length of the sliding window,
        # from full length to a single character
        for length in range(len(sf_a), 0, -1):
            # 'i' stores the starting index of the sliding
            # window in Fx
            for i in range(len(sf_a) - length + 1):
                # extract the pattern for matching
                pattern = sf_a[i : i + length]

                # look for the pattern in Fy
                for idx_b, sf_b in enumerate(f_b):
                    # 'j' stores the starting index in Fy; the
                    # Python find() function returns -1 if there
                    # is no match
                    j = sequence_find(sf_b, pattern)
                    if j is not None:
                        # the pattern was found; set 'new_fx' and
                        # 'new_fy' to version of 'fx' and 'fy' with
                        # the patterns removed, update the SSNC and
                        # set 'match' as True, in order to cascade
                        # out of the loops
                        tmp_x = [sf_a[:i], sf_a[i + length :]]
                        tmp_y = [sf_b[:j], sf_b[j + length :]]
                        new_f_a = f_a[:idx_a] + tmp_x + f_a[idx_a + 1 :]
                        new_f_b = f_b[:idx_b] + tmp_y + f_b[idx_b + 1 :]

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
        new_f_a = [sf for sf in new_f_a if sf]
        new_f_b = [sf for sf in new_f_b if sf]

        return new_f_a, new_f_b, ssnc


def mmcwpa(seq_a: Sequence, seq_b: Sequence) -> float:
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

    # Cache original lengths, before any modifications
    # TODO: make sure it works in a function way, without touching `seq_a`
    #       and `seq_b`, so we don't need to cache this
    len_a, len_b = len(seq_a), len(seq_b)

    # Initialize the `f_a` and `f_b` vectors with a single element each,
    # the corresponding sequence; we also initialize the `ssnc` to zero.
    f_a, f_b = [seq_a], [seq_b]
    ssnc: float = 0.0
    while f_a and f_b:
        f_a, f_b, ssnc = _mmcwpa(f_a, f_b, ssnc)

    return 1.0 - ((ssnc / ((len_a + len_b) ** 2.0)) ** 0.5)
