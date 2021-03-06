"""
Module implementing various methods for similarity and distance from edit methods.

Most of the methods are commonly used in string comparison, such as Levenshtein
distance, but in this module we need to make sure we can operate on arbitrary
iterable data structures.
"""

# Import Python standard libraries
from typing import Callable, Hashable, List, Sequence, Tuple
import difflib
import logging

# Import 3rd-party libraries
import textdistance

# Import local modules
from .common import sequence_find, _nwise, _indices, _wagner_fischer

# Methods for sequence similarity
# -------------------------------
#
# Methods for "similarity" are those that do not guarantee all of the three
# distance properties of positivity, symmetry, identity-discerning, and
# triangle inequality


def birnbaum_simil(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], normal: bool = False
) -> float:
    """
    Compute the Birnbaum similarity score with the fast method.

    This implementation uses the experimental method we developed following
    the description in Birnbaum (2003), which is much faster and less
    memory-intensive than the one implemented in the `birnbaum_simil()`
    function. While in most cases the results are comparable, particularly
    after scaling/normalization, and while the ones
    provided by this method might be considered more adequate due to their
    handling of duplicate information, the values are *not* identical.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.birnbaum_simil("abc", "bcde")
        3.0

    References
    ***********

    Birnbaum, David J. (2003). "Computer-Assisted Analysis and
    Study of the Structure of Mixed-Content Miscellanies". Scripta &
    Scripta 1:15-64.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Whether to normalize the similarity score in range
        [0..1] using sequence lengths.
    :return: The similarity score between the two sequences. The higher
        the similarity score, the more similar the two sequences are;
        a similarity score of zero is the theoretical maximum difference
        between two sequences.
    """

    # If the sequences are equal, we can just compute the score from length
    if seq_x == seq_y:
        length = len(seq_x)
        return (length * (length + 1)) // 2

    # Make sure `seq_x` is shorter or equal in length to `seq_y`
    if len(seq_x) < len(seq_y):
        seq_x, seq_y = seq_y, seq_x

    # Get the opcodes from `SequenceMatcher` and drop the last one (a dummy)
    sm = difflib.SequenceMatcher(None, seq_x, seq_y)
    blocks = sm.get_matching_blocks()[:-1]

    # Sum sizes and return
    sizes = [match.size for match in blocks]
    similarity = sum([(v * (v + 1)) // 2 for v in sizes])

    if normal:
        return float(similarity) / max(
            [
                birnbaum_simil(seq_x, seq_x),
                birnbaum_simil(seq_y, seq_y),
            ]
        )

    return float(similarity)


def fragile_ends_simil(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], normal: bool = False
) -> float:
    """
    Compute the "fragile ends" similarity between two sequences.

    The "fragile ends" similarity is defined as equal to the Levenshtein one,
    but with deletions in the initial or final 10% of the positions being
    discounted.

    This function will use the standard Wagner-Fischer algorithm with the
    default costs provided by the internal `_levdamerau_costs()`
    function. This similarity measure is not used directly in the paper and was a
    proof-of-concept while working toward the "Stemmatological distance".

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.fragile_ends_simil("abc", "bcde")
        3.0

    References
    ***********

    Göransson, Elisabet; Maurits, Luke; Dahlman, Britt; Sarkisian, Karine Å.;
    Rubenson, Samuel; Dunn, Michael. "Improved distance measures for 'mixed-content
    miscellania' (in prep.).

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Whether to normalize the similarity score in range
        [0..1] using sequence lengths.
    :return: The computed "fragile ends" similarity.
    """

    d = _fragile_ends_initial_matrix(seq_x, seq_y)
    similarity = _wagner_fischer(seq_x, seq_y, _fragile_ends_costs, d)

    if normal:
        seq_x_list = [elem for elem in seq_x]
        seq_y_list = [elem for elem in seq_y]

        return similarity / max(
            [
                similarity,
                fragile_ends_simil(seq_x_list + seq_y_list, seq_x_list),
                fragile_ends_simil(seq_x_list + seq_y_list, seq_y_list),
            ]
        )

    return similarity


# TODO: frag type and description
# TODO: should `frag` be in range 0..1, as a float?
def stemmatological_simil(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    frag_start: float = 10.0,
    frag_end: float = 10.0,
    max_del_len: int = 5,
    normal: bool = False,
) -> float:
    """
    Compute the "stemmatological" similarity between two sequences.

    This function will use the standard Wagner-Fischer algorithm with the
    default costs provided by the internal `_levdamerau_costs()`
    function. This similarity measure is essentially a combination of the
    "fragile ends" and "bulk delete" methods, with the first one
    generalised a little to allow specifying the size of both fragile regions.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.stemmatological_simil("abc", "bcde")
        3.0

    References
    ***********

    Göransson, Elisabet; Maurits, Luke; Dahlman, Britt; Sarkisian, Karine Å.;
    Rubenson, Samuel; Dunn, Michael. "Improved distance measures for 'mixed-content
    miscellania' (in prep.).

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param max_del_len: The maximum length of deletion block.
    :param frag_start:
    :param frag_end:
    .param normal: Whether to normalize the similarity score in range
        [0..1] using sequence lengths.
    :return: The computed "stemmatological" similarity.
    """

    d = _stemmatological_initial_matrix(seq_x, seq_y, max_del_len, frag_start, frag_end)
    _stemmatology_costs = _stemmatological_costs_factory(
        max_del_len, frag_start, frag_end
    )

    dist = _wagner_fischer(seq_x, seq_y, _stemmatology_costs, d)

    if normal:
        return dist / max([len(seq_x), len(seq_y)])

    return dist


# Methods for sequence distance
# -----------------------------
#
# Methods for "similarity" are those that guarantee all of the three
# distance properties of positivity, symmetry, identity-discerning, and
# triangle inequality


def levenshtein_dist(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], normal: bool = False
) -> float:
    """
    Compute the Levenshtein distance between two sequences.

    This function will use the standard Wagner-Fischer algorithm with the
    default costs provided by the internal `_levenshtein_costs()`
    function.

    .. code-block:: python

        >>> seqsim.edit.levenshtein_dist("abc", "bcde")
        3.0

    References
    ***********

    Levenshtein, Vladimir I. (February 1966). "Binary codes capable of correcting deletions, insertions,
    and reversals". Soviet Physics Doklady. 10 (8): 707–710

    Wagner, Robert A., and Michael J. Fischer. "The string-to-string correction problem." Journal of the ACM
    (JACM) 21.1 (1974): 168-173.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Whether to normalize the similarity score in range
        [0..1] using sequence lengths.
    :return: The computed Levenshtein distance.
    """

    dist = _wagner_fischer(seq_x, seq_y, _levenshtein_costs)

    if normal:
        return dist / max([len(seq_x), len(seq_y)])

    return dist


def levdamerau_dist(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], normal: bool = False
) -> float:
    """
    Compute the Damerau-Levenshtein distance between two sequences.

    This function will use the standard Wagner-Fischer algorithm with the
    default costs provided by the internal `_levdamerau_costs()`
    function.

    .. code-block:: python

        >>> seqsim.edit.levdamerau_dist("abc", "bcde")
        3.0

    References
    ***********

    Damerau, Fred J. (March 1964), "A technique for computer detection and correction of spelling errors",
    Communications of the ACM, 7 (3): 171–176, doi:10.1145/363958.363994,

    Levenshtein, Vladimir I. (February 1966). "Binary codes capable of correcting deletions, insertions,
    and reversals". Soviet Physics Doklady. 10 (8): 707–710

    Wagner, Robert A., and Michael J. Fischer. "The string-to-string correction problem." Journal of the ACM
    (JACM) 21.1 (1974): 168-173.


    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Whether to normalize the similarity score in range
        [0..1] using sequence lengths.
    :return: The computed Levenshtein distance.
    """

    dist = _wagner_fischer(seq_x, seq_y, _levdamerau_costs)

    if normal:
        return dist / max([len(seq_x), len(seq_y)])

    return dist


# TODO: compute `max_del_len`, if not provided, from sequence length?
def bulk_delete_dist(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    max_del_len: int = 5,
    normal: bool = False,
) -> float:
    """
    Compute the "bulk delete" distance between two sequences.

    This function will use the standard Wagner-Fischer algorithm with the
    default costs provided by the internal `_levdamerau_costs()`
    function. This distance measure is not used directly in the paper and was a
    proof-of-concept while working toward the "Stemmatological distance".

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.bulk_delete_dist("abc", "bcde")
        3.0

    References
    ***********

    Göransson, Elisabet; Maurits, Luke; Dahlman, Britt; Sarkisian, Karine Å.;
    Rubenson, Samuel; Dunn, Michael. "Improved distance measures for 'mixed-content
    miscellania' (in prep.).

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param max_del_len: The maximum length of deletion block.
    :param normal: Whether to normalize the similarity score in range
        [0..1] using sequence lengths.
    :return: The computed "bulk delete" distance.
    """

    d = _bulk_delete_initial_matrix(seq_x, seq_y, max_del_len)
    _bulk_delete_costs = _bulk_delete_costs_factory(max_del_len)
    dist = _wagner_fischer(seq_x, seq_y, _bulk_delete_costs, d)

    if normal:
        return dist / max([len(seq_x), len(seq_y)])

    return dist


def jaro_dist(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], normal: bool = False
) -> float:
    """
    Computes the Jaro distance between two sequences.

    This function returns the value from the implementation provided by
    the `textdistance` library.

    The function accepts the `normal` parameter to have calls equivalent to those
    of other methods, but it is redundant as the Jaccard distance is already
    in range [0..1].

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.jaro_dist("abc", "bcde")
        0.2777777777777778

    References
    ***********

    Jaro, M. A. (1989). "Advances in record linkage methodology as applied to the 1985 census of Tampa Florida".
    Journal of the American Statistical Association. 84 (406): 414–20. doi:10.1080/01621459.1989.10478785.

    Jaro, M. A. (1995). "Probabilistic linkage of large public health data file". Statistics in Medicine. 14 (5–7):
    491–8. doi:10.1002/sim.4780140510. PMID 7792443.

    Winkler, W. E. (1990). "String Comparator Metrics and Enhanced Decision Rules in the Fellegi-Sunter Model of
    Record Linkage". Proceedings of the Section on Survey Research Methods. American Statistical Association: 354–359.

    Winkler, W. E. (2006). "Overview of Record Linkage and Current Research Directions" (PDF). Research
    Report Series, RRS.

    :param seq_x: The first sequence of elements to be compared.
    :param seq_y: The second sequence of elements to be compared.
    :param normal: Dummy parameter, see comment above.
    :return: The Jaro distance between the two sequences.
    """

    dist = textdistance.JaroWinkler(winklerize=False, external=False)(seq_x, seq_y)

    if normal:
        logging.warning(
            "Jaro distance is always in [0..1] range, no need for `normal` parameter."
        )

    return 1.0 - dist


def jaro_winkler_dist(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], normal: bool = False
) -> float:
    """
    Computes the Jaro-Winkler distance between two sequences.

    This function returns the value from the implementation provided by
    the `textdistance` library.

    The function accepts the `normal` parameter to have calls equivalent to those
    of other methods, but it is redundant as the Jaccard distance is already
    in range [0..1].

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.jaro_winkler_dist("abc", "bcde")
        0.2777777777777778

    References
    ***********

    Jaro, M. A. (1989). "Advances in record linkage methodology as applied to the 1985 census of Tampa Florida".
    Journal of the American Statistical Association. 84 (406): 414–20. doi:10.1080/01621459.1989.10478785.

    Jaro, M. A. (1995). "Probabilistic linkage of large public health data file". Statistics in Medicine. 14 (5–7):
    491–8. doi:10.1002/sim.4780140510. PMID 7792443.

    Winkler, W. E. (1990). "String Comparator Metrics and Enhanced Decision Rules in the Fellegi-Sunter Model of
    Record Linkage". Proceedings of the Section on Survey Research Methods. American Statistical Association: 354–359.

    Winkler, W. E. (2006). "Overview of Record Linkage and Current Research Directions" (PDF). Research
    Report Series, RRS.

    :param seq_x: The first sequence of elements to be compared.
    :param seq_y: The second sequence of elements to be compared.
    :param normal: Dummy parameter, see comment above.
    :return: The Jaro-Winkler distance between the two sequences.
    """

    dist = textdistance.JaroWinkler(winklerize=True, external=False)(seq_x, seq_y)

    if normal:
        logging.warning(
            "Jaro-Winkler distance is always in [0..1] range, no need for `normal` parameter."
        )

    return 1.0 - dist


def mmcwpa_dist(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], normal: bool = False
) -> float:
    """
    Computes an MMCWPA distance between two sequences.

    MMCWPA is the Modifier Moving Contracting Window Pattern Algorithm, modified by
    Tiago Tresoldi from a method published by Yang et al. (2001). In order to simplify
    the logic, the function uses the auxiliary internal function `_mmcwpa()`.

    The function accepts the `normal` parameter to have calls equivalent to those
    of other methods, but it is redundant as the Jaccard distance is already
    in range [0..1].

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.mmcwpa_dist("abc", "bcde")
        0.4285714285714286

    References
    ***********

    Tresoldi, Tiago. "Newer method of string comparison: the Modified Moving Contracting Window Pattern Algorithm."
    arXiv preprint arXiv:1605.01079 (2016).

    Yang, Q. X.; Yuan, Sung S.; Chun, Lu; Zhao, Li; Peng Sun. "Faster Algorithm of String
    Comparison", eprint arXiv:cs/0112022, December 2001.

    :param seq_x: The first sequence of elements to be compared.
    :param seq_y: The second sequence of elements to be compared.
    :param normal: Dummy parameter, see comment above.
    :return: The MMCWPA distance between the two sequences.
    """

    # Cache original lengths, before any modifications
    len_a, len_b = len(seq_x), len(seq_y)

    # Initialize the `f_x` and `f_y` vectors with a single element each,
    # the corresponding sequence; we also initialize the `ssnc` to zero.
    f_x, f_y = [seq_x], [seq_y]
    ssnc: float = 0.0
    while f_x and f_y:
        f_x, f_y, ssnc = _mmcwpa(f_x, f_y, ssnc)

    if normal:
        logging.warning(
            "MMCWPA distance is always in [0..1] range, no need for `normal` parameter."
        )

    return 1.0 - ((ssnc / ((len_a + len_b) ** 2.0)) ** 0.5)


def birnbaum_dist(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], normal: bool = False
) -> float:
    """
    Compute the Birnbaum similarity distance with the original method.

    This implementation uses the original method we developed following
    the description in Birnbaum (2003). See comments for `birnbaum_simil()`.

    The function accepts the `normal` parameter to have calls equivalent to those
    of other methods, but it is redundant as the Jaccard distance is already
    in range [0..1].

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.birnbaum_dist("abc", "bcde")
        0.5

    References
    ***********

    Birnbaum, David J. (2003). "Computer-Assisted Analysis and
    Study of the Structure of Mixed-Content Miscellanies". Scripta &
    Scripta 1:15-64.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Dummy parameter, see comment above.
    :return: The distance between the two sequences. A distance of 0.0 indicates
        identical sequences, and a distance of 1.0 indicates the maximum
        theoretical distance between two sequences.
    """

    # Make sure `seq_y` is longer or equal in length to `seq_x`
    if len(seq_y) > len(seq_x):
        seq_x, seq_y = seq_y, seq_x

    distance = 1 - (birnbaum_simil(seq_x, seq_y) / birnbaum_simil(seq_y, seq_y))

    # Compute the distance and correct if necessary (as the distance returned
    # by the method can be less than zero when `seq_x` contains a perfect
    # `seq_y` plus doublets matching parts of `seq_y` as well).
    distance = max(0.0, distance)

    if normal:
        logging.warning(
            "Birnbaum distance is always in [0..1] range, no need for `normal` parameter."
        )

    return distance


def fast_birnbaum_dist(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], normal: bool = False
) -> float:
    """
    Compute the Birnbaum similarity distance with the fast method.

    This implementation uses the experimental method we developed following
    the description in Birnbaum (2003), which is much faster and less
    memory-intensive than the one implemented in the `birnbaum()`
    method. While in most cases the results are comparable, and the ones
    provided by this method might be considered more adequate due to their
    handling of duplicate information, the values are *not* identical.

    The function accepts the `normal` parameter to have calls equivalent to those
    of other methods, but it is redundant as the Jaccard distance is already
    in range [0..1].

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.fast_birnbaum_dist("abc", "bcde")
        0.5

    References
    ***********

    Birnbaum, David J. (2003). "Computer-Assisted Analysis and
    Study of the Structure of Mixed-Content Miscellanies". Scripta &
    Scripta 1:15-64.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Dummy parameter, see comment above.
    :return: The distance between the two sequences. A distance of 0.0 indicates
        identical sequences, and a distance of 1.0 indicates the maximum
        theoretical distance between two sequences.
    """

    # Make sure `seq_y` is longer or equal in length to `seq_x`
    if len(seq_y) > len(seq_x):
        seq_x, seq_y = seq_y, seq_x

    # While the `.birnbaum_simil()` already has an optimization for
    # computing the similarity score of an identity, we can do it here as
    # we know that `seq_y` is equal to itself
    denom = (len(seq_y) * (len(seq_y) + 1)) / 2
    distance = 1.0 - (birnbaum_simil(seq_x, seq_y) / denom)

    if normal:
        logging.warning(
            "Birnbaum distance is always in [0..1] range, no need for `normal` parameter."
        )

    return distance


# Supporting internal functions
# -----------------------------
#
# These are internal functions supporting various different computations,
# mostly related to the computation of matrices for edit distance
# methods. Some of these could (should?) be moved to a new module in
# the future for better organization and a shorter file.


def _levenshtein_costs(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    d: List[List[float]],
    i: int,
    j: int,
) -> Tuple[float, float, float]:
    """
    Computes candidate costs for an entry of an edit distance matrix.

    This internal function will compute the candidate costs for an entry
    (i, j) in terms of the Levenshtein distance matrix (seq_a, seq_b),
    each cost corresponding to one of the available edit operations.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param d: The "starting matrix" for the cost computation.
    :param i: The index of `seq_x` to be considered.
    :param j: The index of `seq_y` to be considered.
    :return: A tuple with the costs for deletion, insertion, and substitution.
    """
    substitution_cost = 0 if seq_x[i - 1] == seq_y[j - 1] else 1
    costs = (
        d[i - 1][j] + 1,  # del
        d[i][j - 1] + 1,  # ins
        d[i - 1][j - 1] + substitution_cost,
    )

    return costs


def _levdamerau_costs(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    d: List[List[float]],
    i: int,
    j: int,
) -> Tuple[float, ...]:
    """
    Computes candidate costs for an entry of a Damerau-Levenshtein distance matrix.

    This internal function will compute the candidate costs for an entry
    (i, j) in terms of the Damerau-Levenshtein distance matrix (seq_a, seq_b),
    each cost corresponding to one of the available edit operations.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param d: The "starting matrix" for the cost computation.
    :param i: The index of `seq_x` to be considered.
    :param j: The index of `seq_y` to be considered.
    :return: A tuple with the costs for deletion, insertion, and substitution,
        and the transposition cost if necessary.
    """

    # Start out as per Levenshtein
    substitution_cost = 0 if seq_x[i - 1] == seq_y[j - 1] else 1
    costs = [
        d[i - 1][j] + 1,
        d[i][j - 1] + 1,
        d[i - 1][j - 1] + substitution_cost,
    ]

    # Add the transposition cost
    if (
        i > 1
        and j > 1
        and seq_x[i - 1] == seq_y[j - 2]
        and seq_x[i - 2] == seq_y[j - 1]
    ):
        costs.append(d[i - 2][j - 2] + 1)

    return tuple(costs)


# TODO: return type and description
def _fragile_ends_initial_matrix(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]
) -> List[List[float]]:
    """
    Compute a starting matrix for the Wagner-Fischer algorithm with "fragile ends" costs.

    The generic starting matrix isn't applicable here, as the cost of deleting the entire
    starting hashable element ("manuscript", in the original application) is less
    than then sequence length.

    See: Göransson, Elisabet; Maurits, Luke; Dahlman, Britt; Sarkisian, Karine Å.;
        Rubenson, Samuel; Dunn, Michael. "Improved distance measures for 'mixed-content
        miscellania' (in prep.).

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :return:
    """

    m = len(seq_x)
    n = len(seq_y)
    d: List[List[float]] = [[0 for i in range(0, n + 1)] for j in range(0, m + 1)]
    for i in range(1, m + 1):
        if i <= round(0.1 * m) or i >= round(0.9 * m):
            d[i][0] = d[i - 1][0] + 0.5
        else:
            d[i][0] = d[i - 1][0] + 1

    # Insertions
    for j in range(1, n + 1):
        d[0][j] = j

    return d


def _fragile_ends_costs(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    d: List[List[float]],
    i: int,
    j: int,
) -> Tuple[float, ...]:
    """
    Computes candidate costs for an entry of a "fragile ends" distance matrix.

    This internal function will compute the candidate costs for an entry
    (i, j) in terms of the Levenshtein distance matrix (seq_a, seq_b),
    each cost corresponding to one of the available edit operations.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param d: The "starting matrix" for the cost computation.
    :param i: The index of `seq_x` to be considered.
    :param j: The index of `seq_y` to be considered.
    :return: A tuple with the costs for deletion, insertion, and substitution.
    """

    substitution_cost = 0 if seq_x[i - 1] == seq_y[j - 1] else 1
    costs = [
        d[i][j - 1] + 1,
        d[i - 1][j - 1] + substitution_cost,
    ]

    # Discount deletion near ends
    m = len(seq_x)
    if i <= round(0.1 * m) or i >= round(0.9 * m):
        costs.append(d[i - 1][j] + 0.5)
    else:
        costs.append(d[i - 1][j] + 1)

    return tuple(costs)


# TODO: return type and description
# TODO: allow max_del_len as a float with percentage length?
def _bulk_delete_initial_matrix(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], max_del_len: int
):
    """
    Compute a starting matrix for the Wagner-Fischer algorithm with "bulk delete" costs.

    The generic starting matrix isn't applicable here, as the cost of deleting the entire
    starting hashable element ("manuscript", in the original application) is less
    than then sequence length.

    See: Göransson, Elisabet; Maurits, Luke; Dahlman, Britt; Sarkisian, Karine Å.;
        Rubenson, Samuel; Dunn, Michael. "Improved distance measures for 'mixed-content
        miscellania' (in prep.).

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param max_del_len: The maximum length of deletion block.
    :return:
    """

    m = len(seq_x)
    n = len(seq_y)
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


# TODO: return description
def _bulk_delete_costs_factory(max_del_len: int = 5) -> Callable:
    """
    Define and return a function for computing candidate costs for a "bulk delete" distance matrix.

    :param max_del_len: The maximum length of deletion block.
    :return:
    """

    def _bulk_delete_costs(
        seq_x: Sequence[Hashable],
        seq_y: Sequence[Hashable],
        d: List[List[float]],
        i: int,
        j: int,
    ) -> Tuple[float, ...]:
        """
        Computes candidate costs for an entry of a "bulk delete" distance matrix.

        This internal function will compute the candidate costs for an entry
        (i, j) in terms of the Levenshtein distance matrix (seq_a, seq_b),
        each cost corresponding to one of the available edit operations.

        :param seq_x: The first sequence to be compared.
        :param seq_y: The second sequence to be compared.
        :param d: The "starting matrix" for the cost computation.
        :param i: The index of `seq_x` to be considered.
        :param j: The index of `seq_y` to be considered.
        :return: A tuple with the costs for deletion, insertion, and substitution.
        """

        substitution_cost = 0 if seq_x[i - 1] == seq_y[j - 1] else 1
        costs = [
            d[i][j - 1] + 1,  # ins
            d[i - 1][j - 1] + substitution_cost,
        ]
        for n in range(1, min(max_del_len + 1, i)):
            # Delete consecutive block of n
            costs.append(d[i - n][j] + 1)

        return tuple(costs)

    return _bulk_delete_costs


# TODO: return type and description
# TODO: frag type and description
def _stemmatological_initial_matrix(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    max_del_len: int = 5,
    frag_start: float = 10.0,
    frag_end: float = 10.0,
):
    """
    Compute a starting matrix for the Wagner-Fischer algorithm with "stemmatological" costs.

    The generic starting matrix isn't applicable here, as the cost of deleting the entire
    starting hashable element ("manuscript", in the original application) is less
    than then sequence length.

    See: Göransson, Elisabet; Maurits, Luke; Dahlman, Britt; Sarkisian, Karine Å.;
        Rubenson, Samuel; Dunn, Michael. "Improved distance measures for 'mixed-content
        miscellania' (in prep.).

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param max_del_len: The maximum length of deletion block.
    :param frag_start:
    :param frag_end:
    :return:
    """
    m = len(seq_x)
    n = len(seq_y)
    d = [[0 for i in range(0, n + 1)] for j in range(0, m + 1)]
    lower = round(m * frag_start / 100.0)
    upper = round(m * (100 - frag_end) / 100.0)
    for i in range(1, m + 1):
        if i <= lower or i >= upper:
            d[i][0]: float = d[i - min(i, max_del_len)][0] + 0.5
        else:
            d[i][0]: float = d[i - min(i, max_del_len)][0] + 1.0

    # Insertions
    for j in range(1, n + 1):
        d[0][j] = j

    return d


# TODO: frag type and description
# TODO: return description
def _stemmatological_costs_factory(
    max_del_len: int = 5, frag_start: float = 10.0, frag_end: float = 10.0
) -> Callable:
    """
    Define and return a function for computing candidate costs for a "stemmatological" distance matrix.

    :param max_del_len: The maximum length of deletion block.
    :param frag_start:
    :param frag_end:
    :return:
    """

    def _stemmatological_costs(
        seq_x: Sequence[Hashable],
        seq_y: Sequence[Hashable],
        d: List[List[float]],
        i: int,
        j: int,
    ):
        """
        Computes candidate costs for an entry of a "stemmatological" distance matrix.

        This internal function will compute the candidate costs for an entry
        (i, j) in terms of the Levenshtein distance matrix (seq_a, seq_b),
        each cost corresponding to one of the available edit operations.

        :param seq_x: The first sequence to be compared.
        :param seq_y: The second sequence to be compared.
        :param d: The "starting matrix" for the cost computation.
        :param i: The index of `seq_x` to be considered.
        :param j: The index of `seq_y` to be considered.
        :return: A tuple with the costs for deletion, insertion, and substitution.
        """
        substitution_cost = 0 if seq_x[i - 1] == seq_y[j - 1] else 1
        costs = [
            d[i][j - 1] + 1,
            d[i - 1][j - 1] + substitution_cost,
        ]
        m = len(seq_x)
        lower = round(m * frag_start / 100.0)
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


def _mmcwpa(
    seq_x: List[Sequence[Hashable]], seq_y: List[Sequence[Hashable]], ssnc: float
) -> Tuple[List[Sequence[Hashable]], List[Sequence[Hashable]], float]:
    """
    Internal function for MMCWPA implementation.

    In this implementation of the Modified Moving Contracting Window Pattern Algorithm
    (MMCWPA) to calculate sequence similarity, we return a list of non-overlapping,
    non-contiguous fields Fx, a list of non-overlapping, non-contiguous fields Fy, and
    the SSNC value (the Sum of the Square of the Number of the same characters). This
    function separates the core method of the implementation and makes recursive calls
    easier.

    :param seq_x: A list of sub-sequences, related to the first sequence.
    :param seq_y: A list of sub-sequences, related to the second sequence.
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
    for idx_x, sf_x in enumerate(seq_x):
        # 'length' stores the length of the sliding window,
        # from full length to a single character
        for length in range(len(sf_x), 0, -1):
            # 'i' stores the starting index of the sliding
            # window in Fx
            for i in range(len(sf_x) - length + 1):
                # extract the pattern for matching
                pattern = sf_x[i : i + length]

                # look for the pattern in Fy
                for idx_y, sf_y in enumerate(seq_y):
                    # 'j' stores the starting index in Fy; the
                    # Python find() function returns -1 if there
                    # is no match
                    j = sequence_find(sf_y, pattern)
                    if j is not None and not match:
                        # the pattern was found; set 'new_fx' and
                        # 'new_fy' to version of 'fx' and 'fy' with
                        # the patterns removed, update the SSNC and
                        # set 'match' as True, in order to cascade
                        # out of the loops
                        tmp_x = [sf_x[:i], sf_x[i + length :]]
                        tmp_y = [sf_y[:j], sf_y[j + length :]]
                        new_f_x = seq_x[:idx_x] + tmp_x + seq_x[idx_x + 1 :]
                        new_f_y = seq_y[:idx_y] + tmp_y + seq_y[idx_y + 1 :]

                        ssnc += (2 * length) ** 2

                        match = True

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
