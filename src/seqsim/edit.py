"""
Module implementing methods for sequence distance and similarity based on edits.

Most of the methods are commonly used in string comparison, such as the
Levenshtein distance, but in this module we make sure we can operate on
arbitrary sequences of hashable elements.

Functions follow a naming convention that states the mathematical properties
of each measure:

  * `_dist`: a true distance (metric), satisfying non-negativity, symmetry,
    identity of indiscernibles (`d(x, y) == 0` if and only if `x == y`), and
    the triangle inequality;
  * `_dissim`: a dissimilarity, where `0.0` is the score for identical
    sequences and higher values indicate more different sequences, but for
    which the metric properties are not all guaranteed;
  * `_simil`: a similarity, where higher values indicate more similar
    sequences.

For the `_dist` edit distances, the properties hold for the raw (i.e., not
normalized) values; normalized values are always in range [0..1], but they
do not satisfy the triangle inequality (e.g., for the Levenshtein distance,
`"ba"`, `"bab"`, and `"ab"`).
"""

# Import Python standard libraries
from typing import Hashable, List, Optional, Sequence
import difflib

# Import local modules
from .common import empty_dissim, sequence_find

# Methods based on the Wagner-Fischer algorithm
# ---------------------------------------------


def levenshtein_dist(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], *, normal: bool = False
) -> float:
    """
    Compute the Levenshtein distance between two sequences.

    The distance is the minimum number of single-element insertions,
    deletions, and substitutions needed to transform one sequence into the
    other, computed with the Wagner-Fischer algorithm.

    Example
    ********

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
    :param normal: Whether to normalize the distance in range [0..1] by
        dividing it by the length of the longest sequence.
    :return: The Levenshtein distance.
    """

    len_y = len(seq_y)
    prev = list(range(len_y + 1))
    for i, elem_x in enumerate(seq_x, start=1):
        curr = [i] + [0] * len_y
        for j, elem_y in enumerate(seq_y, start=1):
            curr[j] = min(
                prev[j] + 1,
                curr[j - 1] + 1,
                prev[j - 1] + (elem_x != elem_y),
            )
        prev = curr

    return _normalize(prev[len_y], seq_x, seq_y, normal)


def osa_dissim(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], *, normal: bool = False
) -> float:
    """
    Compute the Optimal String Alignment (OSA) dissimilarity between two sequences.

    OSA, also known as the "restricted Damerau-Levenshtein distance", extends
    the Levenshtein distance with transpositions of adjacent elements, under
    the restriction that no substring is edited more than once. Due to that
    restriction it does not satisfy the triangle inequality: for example,
    `osa("ca", "abc")` is 3, while `osa("ca", "ac")` + `osa("ac", "abc")` is 2.
    For the unrestricted version, which is a true distance, see
    `damerau_dist()`.

    This was the method offered as `levdamerau_dist()` in previous versions.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.osa_dissim("ca", "abc")
        3.0

    References
    ***********

    Boytsov, Leonid (2011). "Indexing methods for approximate dictionary searching:
    Comparative analysis". Journal of Experimental Algorithmics 16: 1.1.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Whether to normalize the dissimilarity in range [0..1] by
        dividing it by the length of the longest sequence.
    :return: The OSA dissimilarity.
    """

    len_x, len_y = len(seq_x), len(seq_y)
    d = [[0] * (len_y + 1) for _ in range(len_x + 1)]
    for i in range(len_x + 1):
        d[i][0] = i
    for j in range(len_y + 1):
        d[0][j] = j

    for i in range(1, len_x + 1):
        for j in range(1, len_y + 1):
            cost = seq_x[i - 1] != seq_y[j - 1]
            d[i][j] = min(d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + cost)
            if (
                i > 1
                and j > 1
                and seq_x[i - 1] == seq_y[j - 2]
                and seq_x[i - 2] == seq_y[j - 1]
            ):
                d[i][j] = min(d[i][j], d[i - 2][j - 2] + 1)

    return _normalize(d[len_x][len_y], seq_x, seq_y, normal)


def damerau_dist(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], *, normal: bool = False
) -> float:
    """
    Compute the (unrestricted) Damerau-Levenshtein distance between two sequences.

    The distance is the minimum number of insertions, deletions,
    substitutions, and transpositions of two adjacent elements needed to
    transform one sequence into the other, computed with the algorithm by
    Lowrance and Wagner (1975). Unlike the restricted version (see
    `osa_dissim()`), this is a true distance.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.damerau_dist("ca", "abc")
        2.0

    References
    ***********

    Damerau, Fred J. (March 1964), "A technique for computer detection and correction of spelling errors",
    Communications of the ACM, 7 (3): 171–176, doi:10.1145/363958.363994,

    Lowrance, Roy; Wagner, Robert A. (1975). "An Extension of the String-to-String Correction Problem".
    Journal of the ACM 22 (2): 177–183.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Whether to normalize the distance in range [0..1] by
        dividing it by the length of the longest sequence.
    :return: The Damerau-Levenshtein distance.
    """

    len_x, len_y = len(seq_x), len(seq_y)
    max_dist = len_x + len_y

    # The matrix has an extra row and column holding `max_dist`, so that
    # transpositions with elements not seen yet resolve to a sentinel
    d = [[0] * (len_y + 2) for _ in range(len_x + 2)]
    d[0][0] = max_dist
    for i in range(len_x + 1):
        d[i + 1][0] = max_dist
        d[i + 1][1] = i
    for j in range(len_y + 1):
        d[0][j + 1] = max_dist
        d[1][j + 1] = j

    # Last row in which each element was seen in `seq_x`
    last_row = {}
    for i in range(1, len_x + 1):
        elem_x = seq_x[i - 1]
        last_match_col = 0
        for j in range(1, len_y + 1):
            elem_y = seq_y[j - 1]
            k = last_row.get(elem_y, 0)
            ell = last_match_col
            if elem_x == elem_y:
                cost = 0
                last_match_col = j
            else:
                cost = 1
            d[i + 1][j + 1] = min(
                d[i][j] + cost,
                d[i + 1][j] + 1,
                d[i][j + 1] + 1,
                d[k][ell] + (i - k - 1) + 1 + (j - ell - 1),
            )
        last_row[elem_x] = i

    return _normalize(d[len_x + 1][len_y + 1], seq_x, seq_y, normal)


def bulk_delete_dist(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    *,
    max_del_len: int = 5,
    normal: bool = False,
) -> float:
    """
    Compute the "bulk delete" distance between two sequences.

    This is an edit distance where a block of up to `max_del_len` consecutive
    elements can be deleted from, or inserted into, either sequence at the
    cost of a single operation; substitutions have a cost of one. As all
    operations are symmetric and have positive costs, this is a true distance.
    This measure is a proof-of-concept developed while working toward the
    "stemmatological" dissimilarity.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.bulk_delete_dist("abcdeXXXXXfghij", "abcdefghij")
        1.0

    References
    ***********

    Göransson, Elisabet; Maurits, Luke; Dahlman, Britt; Sarkisian, Karine Å.;
    Rubenson, Samuel; Dunn, Michael. "Improved distance measures for 'mixed-content
    miscellania' (in prep.).

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param max_del_len: The maximum length of a block deletion or insertion.
        Must be a positive integer; a value of 1 is equivalent to the
        Levenshtein distance.
    :param normal: Whether to normalize the distance in range [0..1] by
        dividing it by the length of the longest sequence.
    :return: The computed "bulk delete" distance.
    """

    _check_max_del_len(max_del_len)
    dist = _block_edit(seq_x, seq_y, max_del_len, 0.0, 0.0)

    return _normalize(dist, seq_x, seq_y, normal)


def fragile_ends_dissim(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    *,
    frag_start: float = 10.0,
    frag_end: float = 10.0,
    normal: bool = False,
) -> float:
    """
    Compute the "fragile ends" dissimilarity between two sequences.

    The "fragile ends" dissimilarity is equal to the Levenshtein distance, but
    with deletions and insertions in the initial or final positions of each
    sequence (by default, 10% of its length) costing half as much, modelling
    sequences (such as manuscripts) whose beginning and end are more likely to
    be lost. The measure is symmetric, but it does not satisfy the triangle
    inequality, as the discount depends on the positions in each sequence.
    It is a proof-of-concept developed while working toward the
    "stemmatological" dissimilarity.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.fragile_ends_dissim("abcdefghij", "bcdefghij")
        0.5

    References
    ***********

    Göransson, Elisabet; Maurits, Luke; Dahlman, Britt; Sarkisian, Karine Å.;
    Rubenson, Samuel; Dunn, Michael. "Improved distance measures for 'mixed-content
    miscellania' (in prep.).

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param frag_start: The percentage (in range [0..100]) of each sequence,
        from its start, considered fragile.
    :param frag_end: The percentage (in range [0..100]) of each sequence,
        from its end, considered fragile.
    :param normal: Whether to normalize the dissimilarity in range [0..1] by
        dividing it by the length of the longest sequence.
    :return: The computed "fragile ends" dissimilarity.
    """

    _check_frag(frag_start, frag_end)
    dist = _block_edit(seq_x, seq_y, 1, frag_start, frag_end)

    return _normalize(dist, seq_x, seq_y, normal)


def stemmatological_dissim(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    *,
    frag_start: float = 10.0,
    frag_end: float = 10.0,
    max_del_len: int = 5,
    normal: bool = False,
) -> float:
    """
    Compute the "stemmatological" dissimilarity between two sequences.

    This dissimilarity combines the "bulk delete" distance and the "fragile
    ends" dissimilarity: blocks of up to `max_del_len` consecutive elements
    can be deleted or inserted as a single operation, and blocks entirely
    within the fragile regions at the start and end of each sequence cost
    half as much. The measure is symmetric, but it does not satisfy the
    triangle inequality.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.stemmatological_dissim("abcdeXXXXXfghij", "abcdefghij")
        1.0

    References
    ***********

    Göransson, Elisabet; Maurits, Luke; Dahlman, Britt; Sarkisian, Karine Å.;
    Rubenson, Samuel; Dunn, Michael. "Improved distance measures for 'mixed-content
    miscellania' (in prep.).

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param frag_start: The percentage (in range [0..100]) of each sequence,
        from its start, considered fragile.
    :param frag_end: The percentage (in range [0..100]) of each sequence,
        from its end, considered fragile.
    :param max_del_len: The maximum length of a block deletion or insertion.
        Must be a positive integer.
    :param normal: Whether to normalize the dissimilarity in range [0..1] by
        dividing it by the length of the longest sequence.
    :return: The computed "stemmatological" dissimilarity.
    """

    _check_max_del_len(max_del_len)
    _check_frag(frag_start, frag_end)
    dist = _block_edit(seq_x, seq_y, max_del_len, frag_start, frag_end)

    return _normalize(dist, seq_x, seq_y, normal)


# Methods based on matching elements and blocks
# ---------------------------------------------


def jaro_dissim(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], *, normal: bool = False
) -> float:
    """
    Computes the Jaro dissimilarity between two sequences.

    The dissimilarity is one minus the Jaro similarity. As the greedy matching
    of elements depends on the order of the arguments, the similarity is
    computed in both orders and the highest one is used, so that the measure
    is symmetric. It does not satisfy the triangle inequality.

    Results are always in range [0..1], so `normal` has no effect.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.jaro_dissim("abc", "bcde")
        0.2777777777777778

    References
    ***********

    Jaro, M. A. (1989). "Advances in record linkage methodology as applied to the 1985 census of Tampa Florida".
    Journal of the American Statistical Association. 84 (406): 414–20. doi:10.1080/01621459.1989.10478785.

    Jaro, M. A. (1995). "Probabilistic linkage of large public health data file". Statistics in Medicine. 14 (5–7):
    491–8. doi:10.1002/sim.4780140510. PMID 7792443.

    :param seq_x: The first sequence of elements to be compared.
    :param seq_y: The second sequence of elements to be compared.
    :param normal: Ignored, as results are always in range [0..1].
    :return: The Jaro dissimilarity between the two sequences.
    """

    empty = empty_dissim(seq_x, seq_y)
    if empty is not None:
        return empty

    simil = max(
        _jaro_winkler_simil(seq_x, seq_y, winklerize=False),
        _jaro_winkler_simil(seq_y, seq_x, winklerize=False),
    )

    return 1.0 - simil


def jaro_winkler_dissim(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], *, normal: bool = False
) -> float:
    """
    Computes the Jaro-Winkler dissimilarity between two sequences.

    The dissimilarity is one minus the Jaro-Winkler similarity, which boosts
    the Jaro similarity of sequences sharing a common prefix of up to four
    elements. As for `jaro_dissim()`, the similarity is computed in both
    orders and the highest one is used, so that the measure is symmetric. It
    does not satisfy the triangle inequality.

    Results are always in range [0..1], so `normal` has no effect.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.jaro_winkler_dissim("abcd", "abce")
        0.11666666666666659

    References
    ***********

    Winkler, W. E. (1990). "String Comparator Metrics and Enhanced Decision Rules in the Fellegi-Sunter Model of
    Record Linkage". Proceedings of the Section on Survey Research Methods. American Statistical Association: 354–359.

    Winkler, W. E. (2006). "Overview of Record Linkage and Current Research Directions" (PDF). Research
    Report Series, RRS.

    :param seq_x: The first sequence of elements to be compared.
    :param seq_y: The second sequence of elements to be compared.
    :param normal: Ignored, as results are always in range [0..1].
    :return: The Jaro-Winkler dissimilarity between the two sequences.
    """

    empty = empty_dissim(seq_x, seq_y)
    if empty is not None:
        return empty

    simil = max(
        _jaro_winkler_simil(seq_x, seq_y, winklerize=True),
        _jaro_winkler_simil(seq_y, seq_x, winklerize=True),
    )

    return 1.0 - simil


def mmcwpa_dissim(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], *, normal: bool = False
) -> float:
    """
    Computes the MMCWPA dissimilarity between two sequences.

    MMCWPA is the Modified Moving Contracting Window Pattern Algorithm,
    modified by Tiago Tresoldi from a method published by Yang et al. (2001).
    It repeatedly removes the longest common sub-sequence (found with a
    contracting window) from both sequences, accumulating the SSNC ("sum of
    the squares of the number of the same characters"), and returns
    `1 - sqrt(SSNC / (len(x) + len(y)) ** 2)`. As the greedy search depends on
    the order of the arguments, it is performed in both orders and the
    highest SSNC is used, so that the measure is symmetric. It does not
    satisfy the triangle inequality.

    Results are always in range [0..1], so `normal` has no effect.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.mmcwpa_dissim("abc", "bcde")
        0.4285714285714286

    References
    ***********

    Tresoldi, Tiago. "Newer method of string comparison: the Modified Moving Contracting Window Pattern Algorithm."
    arXiv preprint arXiv:1605.01079 (2016).

    Yang, Q. X.; Yuan, Sung S.; Chun, Lu; Zhao, Li; Peng Sun. "Faster Algorithm of String
    Comparison", eprint arXiv:cs/0112022, December 2001.

    :param seq_x: The first sequence of elements to be compared.
    :param seq_y: The second sequence of elements to be compared.
    :param normal: Ignored, as results are always in range [0..1].
    :return: The MMCWPA dissimilarity between the two sequences.
    """

    empty = empty_dissim(seq_x, seq_y)
    if empty is not None:
        return empty

    ssnc = max(_mmcwpa_ssnc(seq_x, seq_y), _mmcwpa_ssnc(seq_y, seq_x))

    return 1.0 - ((ssnc / ((len(seq_x) + len(seq_y)) ** 2.0)) ** 0.5)


def birnbaum_simil(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], *, normal: bool = False
) -> float:
    """
    Compute the Birnbaum similarity score between two sequences.

    This implementation follows the description in Birnbaum (2003): the
    sequences are aligned into matching blocks (with the Ratcliff-Obershelp
    algorithm of Python's `difflib`), and each block of size `n` contributes
    `n * (n + 1) / 2` to the score, so that longer blocks weigh more. As the
    alignment depends on the order of the arguments, it is computed in both
    orders and the highest score is used, so that the measure is symmetric.

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
    :param normal: Whether to normalize the similarity score in range [0..1]
        by dividing it by the score of the longest sequence compared with
        itself. A normalized score of 1.0 indicates identical sequences.
    :return: The similarity score between the two sequences. The higher the
        score, the more similar the two sequences are; a score of zero
        indicates that no element is shared.
    """

    max_len = max(len(seq_x), len(seq_y))
    if max_len == 0:
        return 1.0 if normal else 0.0

    similarity = max(_birnbaum_score(seq_x, seq_y), _birnbaum_score(seq_y, seq_x))

    if normal:
        return similarity / _triangular(max_len)

    return float(similarity)


def birnbaum_dissim(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], *, normal: bool = False
) -> float:
    """
    Compute the Birnbaum dissimilarity between two sequences.

    The dissimilarity is one minus the normalized Birnbaum similarity (see
    `birnbaum_simil()`), so that 0.0 indicates identical sequences and 1.0
    sequences sharing no element. It does not satisfy the triangle
    inequality.

    Results are always in range [0..1], so `normal` has no effect.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.birnbaum_dissim("abc", "bcde")
        0.7

    References
    ***********

    Birnbaum, David J. (2003). "Computer-Assisted Analysis and
    Study of the Structure of Mixed-Content Miscellanies". Scripta &
    Scripta 1:15-64.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Ignored, as results are always in range [0..1].
    :return: The Birnbaum dissimilarity between the two sequences.
    """

    return 1.0 - birnbaum_simil(seq_x, seq_y, normal=True)


# Supporting internal functions
# -----------------------------


def _normalize(
    dist: float, seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], normal: bool
) -> float:
    """
    Returns an edit distance as a float, normalized by the longest length if requested.

    For all edit measures in this module, the length of the longest sequence
    is an upper bound (it is the cost of substituting and inserting or
    deleting elements one by one), so normalized values are in range [0..1].
    Two empty sequences have a normalized distance of 0.0.
    """

    if normal:
        max_len = max(len(seq_x), len(seq_y))
        return dist / max_len if max_len else 0.0

    return float(dist)


def _check_max_del_len(max_del_len: int) -> None:
    """
    Raises a `ValueError` if `max_del_len` is not a positive integer.
    """

    if isinstance(max_del_len, bool) or not isinstance(max_del_len, int):
        raise ValueError(f"`max_del_len` must be an integer, got {max_del_len!r}.")
    if max_del_len < 1:
        raise ValueError(f"`max_del_len` must be at least 1, got {max_del_len}.")


def _check_frag(frag_start: float, frag_end: float) -> None:
    """
    Raises a `ValueError` if the fragile region percentages are out of range.
    """

    for name, value in (("frag_start", frag_start), ("frag_end", frag_end)):
        if not 0.0 <= value <= 100.0:
            raise ValueError(f"`{name}` must be in range [0..100], got {value!r}.")


def _fragile_bounds(length: int, frag_start: float, frag_end: float):
    """
    Returns the bounds of the fragile regions of a sequence.

    The first `lower` positions and all positions after `upper` (with
    positions counted from 1) are fragile.
    """

    lower = round(length * frag_start / 100.0)
    upper = round(length * (100.0 - frag_end) / 100.0)

    return lower, upper


def _block_edit(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    max_block: int,
    frag_start: float,
    frag_end: float,
    fragile_cost: float = 0.5,
) -> float:
    """
    Computes an edit distance with block operations and fragile ends.

    Blocks of up to `max_block` consecutive elements can be deleted from
    `seq_x` or inserted from `seq_y` at a cost of 1.0, or `fragile_cost` if
    the block lies entirely within a fragile region of its sequence (at the
    start or end, as percentages of its length). Substitutions cost 1.0. As
    deletions from `seq_x` and insertions from `seq_y` follow the same rules,
    the result is symmetric.
    """

    len_x, len_y = len(seq_x), len(seq_y)
    lower_x, upper_x = _fragile_bounds(len_x, frag_start, frag_end)
    lower_y, upper_y = _fragile_bounds(len_y, frag_start, frag_end)

    def block_cost(start: int, end: int, lower: int, upper: int) -> float:
        # Block of positions `start + 1` to `end` (1-based, inclusive)
        return fragile_cost if end <= lower or start >= upper else 1.0

    d: List[List[float]] = [[0.0] * (len_y + 1) for _ in range(len_x + 1)]
    for i in range(len_x + 1):
        for j in range(len_y + 1):
            if i == 0 and j == 0:
                continue

            candidates = []
            if i and j:
                candidates.append(d[i - 1][j - 1] + (seq_x[i - 1] != seq_y[j - 1]))
            for size in range(1, min(max_block, i) + 1):
                candidates.append(
                    d[i - size][j] + block_cost(i - size, i, lower_x, upper_x)
                )
            for size in range(1, min(max_block, j) + 1):
                candidates.append(
                    d[i][j - size] + block_cost(j - size, j, lower_y, upper_y)
                )

            d[i][j] = min(candidates)

    return d[len_x][len_y]


def _jaro_winkler_simil(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    winklerize: bool,
    prefix_weight: float = 0.1,
) -> float:
    """
    Computes the Jaro (or Jaro-Winkler) similarity between two sequences.

    This follows the implementation in the `textdistance` library (version
    4.5), from which it was ported: elements of `seq_x` are matched to the
    first unmatched equal element of `seq_y` within the search window, and the
    Winkler boost is only applied when the Jaro similarity is above 0.7.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param winklerize: Whether to apply the Winkler prefix boost.
    :param prefix_weight: The scaling factor for the Winkler prefix boost.
    :return: The similarity, in range [0..1].
    """

    len_x, len_y = len(seq_x), len(seq_y)
    if tuple(seq_x) == tuple(seq_y):
        return 1.0
    if not len_x or not len_y:
        return 0.0

    search_range = max(max(len_x, len_y) // 2 - 1, 0)

    # Flag matching elements within the search range
    flags_x = [False] * len_x
    flags_y = [False] * len_y
    common = 0
    for i, elem_x in enumerate(seq_x):
        low = max(0, i - search_range)
        high = min(i + search_range, len_y - 1)
        for j in range(low, high + 1):
            if not flags_y[j] and seq_y[j] == elem_x:
                flags_x[i] = flags_y[j] = True
                common += 1
                break

    if not common:
        return 0.0

    # Count transpositions
    k = transpositions = 0
    for i, flag_x in enumerate(flags_x):
        if flag_x:
            for j in range(k, len_y):
                if flags_y[j]:
                    k = j + 1
                    break
            if seq_x[i] != seq_y[j]:
                transpositions += 1
    transpositions //= 2

    simil = (common / len_x + common / len_y + (common - transpositions) / common) / 3

    # Apply the Winkler boost for a common prefix of up to four elements
    if winklerize and simil > 0.7:
        prefix = 0
        while prefix < min(len_x, len_y, 4) and seq_x[prefix] == seq_y[prefix]:
            prefix += 1
        simil += prefix * prefix_weight * (1.0 - simil)

    return simil


def _triangular(value: int) -> int:
    """
    Returns the triangular number of `value`.
    """

    return (value * (value + 1)) // 2


def _birnbaum_score(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> int:
    """
    Computes the (directional) Birnbaum similarity score.

    The automatic junk heuristic of `SequenceMatcher` is disabled, as it
    would otherwise change the results for sequences with 200 or more
    elements.
    """

    matcher = difflib.SequenceMatcher(None, seq_x, seq_y, autojunk=False)

    # The last matching block is a dummy of size zero
    return sum(_triangular(block.size) for block in matcher.get_matching_blocks())


def _mmcwpa_ssnc(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
    """
    Computes the (directional) SSNC of the MMCWPA method.
    """

    f_x: List[Sequence[Hashable]] = [seq_x]
    f_y: List[Sequence[Hashable]] = [seq_y]
    ssnc = 0.0
    while f_x and f_y:
        f_x, f_y, ssnc = _mmcwpa(f_x, f_y, ssnc)

    return ssnc


def _mmcwpa(seq_x, seq_y, ssnc):
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
             from the second sequence, and the third element is the updated SSNC. If
             no match is found, both lists are empty.
    """

    # Search patterns in all subfields of Fx, in order, using a window
    # contracting from the full length of the subfield to a single element;
    # the first match found is removed from both Fx and Fy
    for idx_x, sf_x in enumerate(seq_x):
        for length in range(len(sf_x), 0, -1):
            for i in range(len(sf_x) - length + 1):
                pattern = sf_x[i : i + length]
                for idx_y, sf_y in enumerate(seq_y):
                    j: Optional[int] = sequence_find(sf_y, pattern)
                    if j is not None:
                        new_f_x = (
                            seq_x[:idx_x]
                            + [sf_x[:i], sf_x[i + length :]]
                            + seq_x[idx_x + 1 :]
                        )
                        new_f_y = (
                            seq_y[:idx_y]
                            + [sf_y[:j], sf_y[j + length :]]
                            + seq_y[idx_y + 1 :]
                        )

                        # Remove any empty subfields due to pattern removal
                        new_f_x = [sf for sf in new_f_x if len(sf)]
                        new_f_y = [sf for sf in new_f_y if len(sf)]

                        return new_f_x, new_f_y, ssnc + (2 * length) ** 2

    return [], [], ssnc
