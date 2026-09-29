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
from typing import Hashable, List, Sequence
import difflib

# Import local modules
from ._measure import Scored, measure
from .common import equivalent_string, lcs_length


class _Boundary:
    """
    Sentinel for sequence boundaries.
    """

    def __init__(self, name: str):
        self.name = name

    def __repr__(self) -> str:
        return self.name


_BLOCK_START = _Boundary("START")
_BLOCK_END = _Boundary("END")

# Validation and normalization helpers
# -------------------------------------


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


def _check_block_options(
    max_del_len: int = 1, frag_start: float = 0.0, frag_end: float = 0.0
) -> None:
    """
    Validates the options of the block edit measures.
    """

    _check_max_del_len(max_del_len)
    _check_frag(frag_start, frag_end)


def _check_min_match(min_match: int) -> None:
    """
    Raises a `ValueError` if `min_match` is not a positive integer.
    """

    if isinstance(min_match, bool) or not isinstance(min_match, int) or min_match < 1:
        raise ValueError(f"`min_match` must be a positive integer, got {min_match!r}.")


def _gld(dist: float, seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
    """
    Returns the normalization of Yujian and Bo (2007) of an edit distance.

    With unit insertion and deletion costs, the normalized value is
    `2 * d / (len(x) + len(y) + d)`.
    """

    denominator = len(seq_x) + len(seq_y) + dist

    return 2.0 * dist / denominator if denominator else 0.0


# Methods based on the Wagner-Fischer algorithm
# ---------------------------------------------


@measure(key="levenshtein", kind="dist", triangle="yes", bound="max_len")
def levenshtein_dist(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
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

    return prev[len_y]


@measure(key="levenshtein_gld", kind="dist", triangle="yes")
def levenshtein_gld_dist(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
    """
    Compute the normalized Levenshtein distance of Yujian and Bo (2007).

    The normalized distance is `2 * d / (len(x) + len(y) + d)`, where `d` is
    the Levenshtein distance. Unlike dividing `d` by the length of the longest
    sequence, this normalization is proven to satisfy the triangle
    inequality, so the result is a true distance in range [0..1]. It is 1.0
    only when one of the sequences is empty.

    Results are always in range [0..1], so `normal` has no effect.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.levenshtein_gld_dist("kitten", "sitting")
        0.375

    References
    ***********

    Yujian, Li; Bo, Liu (2007). "A Normalized Levenshtein Distance Metric". IEEE
    Transactions on Pattern Analysis and Machine Intelligence 29 (6): 1091–1095.
    doi:10.1109/TPAMI.2007.1078

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :return: The normalized Levenshtein distance.
    """

    return _gld(levenshtein_dist(seq_x, seq_y), seq_x, seq_y)


@measure(key="levenshtein_ned", kind="dist", triangle="yes")
def levenshtein_ned_dist(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
    """
    Compute the normalized edit distance of Marzal and Vidal (1993).

    The normalized edit distance is the minimum, over all edit paths
    transforming one sequence into the other, of the number of edits
    (insertions, deletions, and substitutions, each with cost one) divided by
    the length of the path (the number of operations, including the matches
    of equal elements). This is not the same as dividing the Levenshtein
    distance by a length, as a longer path with more matches can have a lower
    ratio. With unit costs it is proven to be a true distance (Fisman et al.,
    2022) in range [0..1].

    It is computed in `O(len(x) * len(y) * (len(x) + len(y)))` time, so it is
    considerably slower than the other edit distances for long sequences.
    Results are always in range [0..1], so `normal` has no effect.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.levenshtein_ned_dist("kitten", "sitting")
        0.42857142857142855

    References
    ***********

    Marzal, Andrés; Vidal, Enrique (1993). "Computation of normalized edit distance
    and applications". IEEE Transactions on Pattern Analysis and Machine
    Intelligence 15 (9): 926–932. doi:10.1109/34.232078

    Fisman, Dana; Grogin, Joshua; Margalit, Oded; Weiss, Gera (2022). "The Normalized
    Edit Distance with Uniform Operation Costs is a Metric". 33rd Annual Symposium on
    Combinatorial Pattern Matching (CPM 2022), LIPIcs 223: 17:1–17:17.
    doi:10.4230/LIPIcs.CPM.2022.17

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :return: The normalized edit distance.
    """

    len_x, len_y = len(seq_x), len(seq_y)
    if not len_x and not len_y:
        return 0.0

    # `layer[i][j]` holds the minimum number of edits of a path from (0, 0) to
    # (i, j) with exactly `steps` operations
    inf = float("inf")
    layer = [[inf] * (len_y + 1) for _ in range(len_x + 1)]
    layer[0][0] = 0
    best = inf
    for steps in range(1, len_x + len_y + 1):
        new = [[inf] * (len_y + 1) for _ in range(len_x + 1)]
        # A path of `steps` operations can only reach cells with
        # max(i, j) <= steps <= i + j
        for i in range(min(steps, len_x) + 1):
            for j in range(max(0, steps - i), min(steps, len_y) + 1):
                candidates = []
                if i and j:
                    candidates.append(
                        layer[i - 1][j - 1] + (seq_x[i - 1] != seq_y[j - 1])
                    )
                if i:
                    candidates.append(layer[i - 1][j] + 1)
                if j:
                    candidates.append(layer[i][j - 1] + 1)
                new[i][j] = min(candidates)
        layer = new
        if layer[len_x][len_y] < inf:
            best = min(best, layer[len_x][len_y] / steps)

    return best


@measure(
    key="osa",
    kind="dissim",
    triangle="no",
    triangle_example=("ca", "ac", "abc"),
    bound="max_len",
)
def osa_dissim(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
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

    return d[len_x][len_y]


@measure(key="damerau", kind="dist", triangle="yes", bound="max_len")
def damerau_dist(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
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

    return d[len_x + 1][len_y + 1]


@measure(key="indel", kind="dist", triangle="yes", bound="sum_len")
def indel_dist(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
    """
    Compute the insertion-deletion (indel) distance between two sequences.

    The distance is the minimum number of single-element insertions and
    deletions needed to transform one sequence into the other, that is,
    `len(x) + len(y) - 2 * LCS(x, y)`, where `LCS` is the length of the
    longest common (not necessarily contiguous) subsequence. It is an edit
    distance with symmetric, unit-cost operations, and thus a true distance.
    It is a natural measure for lists of contents, where each element is
    either shared or not, in order.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.indel_dist("kitten", "sitting")
        5.0

    References
    ***********

    Needleman, Saul B.; Wunsch, Christian D. (1970). "A general method applicable to
    the search for similarities in the amino acid sequence of two proteins". Journal
    of Molecular Biology 48 (3): 443–53.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :return: The indel distance.
    """

    return len(seq_x) + len(seq_y) - 2 * lcs_length(seq_x, seq_y)


@measure(key="lcs", kind="dist", triangle="yes")
def lcs_dist(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
    """
    Compute the normalized longest common subsequence (LCS) distance.

    The distance is `1 - LCS(x, y) / max(len(x), len(y))`, where `LCS` is the
    length of the longest common (not necessarily contiguous) subsequence. It
    is the proportion of the longest sequence not covered by the common
    subsequence, and it is a true distance in range [0..1] (Bakkelund, 2009).

    Results are always in range [0..1], so `normal` has no effect.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.lcs_dist("kitten", "sitting")
        0.4285714285714286

    References
    ***********

    Bakkelund, Daniel (2009). "An LCS-based string metric". Technical report,
    University of Oslo.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :return: The LCS distance.
    """

    max_len = max(len(seq_x), len(seq_y))
    if not max_len:
        return 0.0

    return 1.0 - lcs_length(seq_x, seq_y) / max_len


@measure(key="damerau_gld", kind="dist", triangle="yes")
def damerau_gld_dist(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
    """
    Compute the normalized Damerau-Levenshtein distance of Yujian and Bo (2007).

    This is the normalization of `levenshtein_gld_dist()`,
    `2 * d / (len(x) + len(y) + d)`, applied to the (unrestricted)
    Damerau-Levenshtein distance `d` (see `damerau_dist()`). The formula is
    the Steinhaus transform of `d` with respect to the empty sequence, which
    preserves the metric properties for any metric `d` for which the
    distance of a sequence to the empty sequence is its length; the result
    is thus a true distance in range [0..1].

    Results are always in range [0..1], so `normal` has no effect.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.damerau_gld_dist("ca", "abc")
        0.5714285714285714

    References
    ***********

    Yujian, Li; Bo, Liu (2007). "A Normalized Levenshtein Distance Metric". IEEE
    Transactions on Pattern Analysis and Machine Intelligence 29 (6): 1091–1095.
    doi:10.1109/TPAMI.2007.1078

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :return: The normalized Damerau-Levenshtein distance.
    """

    return _gld(damerau_dist(seq_x, seq_y), seq_x, seq_y)


@measure(key="indel_gld", kind="dist", triangle="yes")
def indel_gld_dist(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
    """
    Compute the normalized indel distance of Yujian and Bo (2007).

    This is the normalization of `levenshtein_gld_dist()`,
    `2 * d / (len(x) + len(y) + d)`, applied to the insertion-deletion
    distance `d` (see `indel_dist()`). As for `damerau_gld_dist()`, the
    result is a true distance in range [0..1].

    Results are always in range [0..1], so `normal` has no effect.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.indel_gld_dist("kitten", "sitting")
        0.5555555555555556

    References
    ***********

    Yujian, Li; Bo, Liu (2007). "A Normalized Levenshtein Distance Metric". IEEE
    Transactions on Pattern Analysis and Machine Intelligence 29 (6): 1091–1095.
    doi:10.1109/TPAMI.2007.1078

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :return: The normalized indel distance.
    """

    return _gld(indel_dist(seq_x, seq_y), seq_x, seq_y)


@measure(
    key="bulk_delete",
    kind="dist",
    triangle="yes",
    bound="max_len",
    check=_check_block_options,
)
def bulk_delete_dist(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    *,
    max_del_len: int = 5,
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
    :return: The computed "bulk delete" distance.
    """

    return _block_edit(seq_x, seq_y, max_del_len, 0.0, 0.0)


@measure(
    key="fragile_ends",
    kind="dissim",
    triangle="no",
    triangle_example=("baabb", "bbaababbba", "bbbaabbabbbaa"),
    bound="max_len",
    check=_check_block_options,
)
def fragile_ends_dissim(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    *,
    frag_start: float = 10.0,
    frag_end: float = 10.0,
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
    :return: The computed "fragile ends" dissimilarity.
    """

    return _block_edit(seq_x, seq_y, 1, frag_start, frag_end)


@measure(
    key="stemmatological",
    kind="dissim",
    triangle="no",
    triangle_example=("baba", "babab", "ab"),
    bound="max_len",
    check=_check_block_options,
)
def stemmatological_dissim(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    *,
    frag_start: float = 10.0,
    frag_end: float = 10.0,
    max_del_len: int = 5,
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
    :return: The computed "stemmatological" dissimilarity.
    """

    return _block_edit(seq_x, seq_y, max_del_len, frag_start, frag_end)


# Methods based on matching elements and blocks
# ---------------------------------------------


@measure(
    key="jaro",
    kind="dissim",
    triangle="no",
    triangle_example=("baaa", "cba", "cccc"),
    empty="max",
    symmetrize="min",
)
def jaro_dissim(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
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
    :return: The Jaro dissimilarity between the two sequences.
    """

    return 1.0 - _jaro_winkler_simil(seq_x, seq_y, winklerize=False)


@measure(
    key="jaro_winkler",
    kind="dissim",
    triangle="no",
    triangle_example=("caac", "bbba", "bbb"),
    empty="max",
    symmetrize="min",
)
def jaro_winkler_dissim(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
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
    :return: The Jaro-Winkler dissimilarity between the two sequences.
    """

    return 1.0 - _jaro_winkler_simil(seq_x, seq_y, winklerize=True)


@measure(
    key="mmcwpa",
    kind="dissim",
    triangle="no",
    triangle_example=("acca", "cabb", "bbb"),
    empty="max",
    symmetrize="min",
)
def mmcwpa_dissim(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
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
    :return: The MMCWPA dissimilarity between the two sequences.
    """

    ssnc = _mmcwpa_ssnc(seq_x, seq_y)

    return 1.0 - ((ssnc / ((len(seq_x) + len(seq_y)) ** 2.0)) ** 0.5)


@measure(
    kind="simil",
    bound="custom",
    symmetrize="max",
    raw_range="0 upwards",
    normal_doc=(
        "Whether to normalize the similarity score in range [0..1] by "
        "dividing it by the score of the longest sequence compared with "
        "itself. A normalized score of 1.0 indicates identical sequences."
    ),
)
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
    :return: The similarity score between the two sequences. The higher the
        score, the more similar the two sequences are; a score of zero
        indicates that no element is shared.
    """

    max_len = max(len(seq_x), len(seq_y))
    if max_len == 0:
        return 1.0 if normal else 0.0

    similarity = _birnbaum_score(seq_x, seq_y)
    if normal:
        return similarity / _triangular(max_len)

    return similarity


@measure(
    key="birnbaum",
    kind="dissim",
    triangle="no",
    triangle_example=("bbbaaabbaaaa", "bbbaaaabaaaa", "bbaaaabaaaa"),
)
def birnbaum_dissim(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
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
    :return: The Birnbaum dissimilarity between the two sequences.
    """

    return 1.0 - birnbaum_simil(seq_x, seq_y, normal=True)


# Methods based on block moves
# ----------------------------


@measure(
    key="block_move",
    kind="dissim",
    triangle="no",
    triangle_example=("a", "aa", "aaaa"),
    bound="scored",
    raw_range="0 to max length + 1",
    symmetrize="max",
    directional_option=True,
    normal_doc=(
        "Whether to normalize the result in range [0..1] by dividing it by "
        "one plus the length of the longest sequence, which is an upper bound."
    ),
)
def block_move_dissim(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
) -> float:
    """
    Computes the block move dissimilarity between two sequences.

    Following Tichy (1984), a sequence `y` can be built from `x` by copying
    blocks (sub-sequences) of `x`, in any order and possibly more than once,
    and adding the elements of `y` not found in `x`. The minimum number of
    pieces (block moves, plus one addition for each maximal run of elements
    not found in `x`) is found greedily, by taking at each
    position of `y` the longest prefix of its remainder found in `x`, which
    Tichy proves to be optimal. Both sequences are wrapped in start and end
    boundaries, so that identical sequences need a single block, and the
    dissimilarity is the number of pieces minus one: the number of "cuts"
    needed to build one sequence from the other. As this is directional, the
    highest number for both directions is returned, unless `directional` is
    set. Different sequences have a positive dissimilarity, but the triangle
    inequality does not always hold.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.block_move_dissim("abcdefgh", "efghabcd")
        3.0
        >>> seqsim.edit.block_move_dissim("abcdefgh", "abcdXefgh")
        2.0

    References
    ***********

    Tichy, Walter F. (1984). "The string-to-string correction problem with block
    moves". ACM Transactions on Computer Systems 2 (4): 309–321.
    doi:10.1145/357401.357404

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param directional: Whether to return the number of cuts needed to build
        `seq_y` from `seq_x` only. Defaults to `False`.
    :return: The block move dissimilarity.
    """

    str_x, str_y = equivalent_string(
        [_BLOCK_START, *seq_x, _BLOCK_END], [_BLOCK_START, *seq_y, _BLOCK_END]
    )
    cuts = _block_cover(str_x, str_y) - 1

    return Scored(cuts, max(len(seq_x), len(seq_y)) + 1)


def _block_cover(source: str, target: str) -> int:
    """
    Returns the minimum number of pieces needed to build `target` from `source`.
    """

    pieces = 0
    pos = 0
    while pos < len(target):
        length = 0
        while pos + length < len(target) and target[pos : pos + length + 1] in source:
            length += 1
        if not length:
            # Elements not found in `source` are added, one maximal run of
            # such elements at a time
            while pos + length < len(target) and target[pos + length] not in source:
                length += 1
        pos += length
        pieces += 1

    return pieces


@measure(
    key="gst",
    kind="dissim",
    identity="no",
    identity_example=("aaab", "abaa"),
    triangle="no",
    triangle_example=("aa", "aab", "ab"),
    empty="max",
    symmetrize="min",
    check=_check_min_match,
)
def gst_dissim(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    *,
    min_match: int = 2,
) -> float:
    """
    Computes the Greedy String Tiling dissimilarity between two sequences.

    Greedy String Tiling (Wise, 1993) covers both sequences with "tiles",
    non-overlapping common sub-sequences of at least `min_match` elements,
    taking the longest available matches first. As tiles can be found in any
    order, it tolerates moved blocks (such as transposed groups of texts),
    unlike edit distances. The similarity is the proportion of elements
    covered by tiles, `2 * coverage / (len(x) + len(y))`, and the
    dissimilarity is one minus it. As ties between matches of the same length
    depend on the order of the arguments, the tiling is computed in both
    orders and the highest coverage is used. Identical sequences always
    have a dissimilarity of zero, even if shorter than `min_match`. The
    measure ignores the order of the tiles, so different sequences can have a dissimilarity of zero (e.g.,
    `"abcd"` and `"cdab"` with `min_match=2`).

    Results are always in range [0..1], so `normal` has no effect.

    Example
    ********

    .. code-block:: python

        >>> seqsim.edit.gst_dissim("abcdefgh", "efghabcd")
        0.0
        >>> seqsim.edit.gst_dissim("abcdefgh", "efghXbcd")
        0.125

    References
    ***********

    Wise, Michael J. (1993). "String Similarity via Greedy String Tiling and Running
    Karp-Rabin Matching". Technical report, Department of Computer Science,
    University of Sydney.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param min_match: The minimum length of a tile. Defaults to 2, so that
        single shared elements out of context are not counted.
    :return: The Greedy String Tiling dissimilarity.
    """

    # Identical sequences shorter than `min_match` cannot be tiled
    if tuple(seq_x) == tuple(seq_y):
        return 0.0

    str_x, str_y = equivalent_string(seq_x, seq_y)
    coverage = _gst_coverage(str_x, str_y, min_match)

    return 1.0 - (2.0 * coverage / (len(str_x) + len(str_y)))


def _gst_coverage(str_x: str, str_y: str, min_match: int) -> int:
    """
    Returns the number of elements of `str_x` covered by Greedy String Tiling.
    """

    marked_x = [False] * len(str_x)
    marked_y = [False] * len(str_y)
    coverage = 0

    while True:
        # Find all maximal matches of the longest length, among unmarked
        # elements
        max_match = min_match
        matches: List = []
        for i in range(len(str_x)):
            if marked_x[i]:
                continue
            for j in range(len(str_y)):
                length = 0
                while (
                    i + length < len(str_x)
                    and j + length < len(str_y)
                    and not marked_x[i + length]
                    and not marked_y[j + length]
                    and str_x[i + length] == str_y[j + length]
                ):
                    length += 1
                if length > max_match:
                    matches = [(i, j, length)]
                    max_match = length
                elif length == max_match:
                    matches.append((i, j, length))

        # Mark the tiles that are not occluded by tiles marked earlier
        for i, j, length in matches:
            if not any(marked_x[i : i + length]) and not any(marked_y[j : j + length]):
                marked_x[i : i + length] = [True] * length
                marked_y[j : j + length] = [True] * length
                coverage += length

        if max_match == min_match:
            return coverage


# Supporting internal functions
# -----------------------------


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

    The sequences are lists of non-overlapping, non-contiguous subfields Fx
    and Fy, initially holding the full sequences. At each step, the
    subfields of Fx are searched in order with a window contracting from the
    full length of the subfield to a single element, and the first pattern
    found in a subfield of Fy (searched in order) is removed from both,
    adding the square of the number of matching elements to the SSNC ("Sum
    of the Square of the Number of the same characters"). The search ends
    when no pattern is found.

    For efficiency, sequences are first mapped to equivalent strings, and the
    longest window with a match in each subfield of Fx is found with a binary
    search (as a match of length `n` implies matches of all shorter lengths),
    which is equivalent to trying all window lengths in decreasing order.
    """

    str_x, str_y = equivalent_string(seq_x, seq_y)

    f_x: List[str] = [str_x]
    f_y: List[str] = [str_y]
    ssnc = 0.0
    while f_x and f_y:
        match = _mmcwpa_find(f_x, f_y)
        if match is None:
            break

        idx_x, i, idx_y, j, length = match
        sf_x, sf_y = f_x[idx_x], f_y[idx_y]
        f_x[idx_x : idx_x + 1] = [sf for sf in (sf_x[:i], sf_x[i + length :]) if sf]
        f_y[idx_y : idx_y + 1] = [sf for sf in (sf_y[:j], sf_y[j + length :]) if sf]
        ssnc += (2 * length) ** 2

    return ssnc


def _mmcwpa_find(f_x: List[str], f_y: List[str]):
    """
    Finds the next match of the MMCWPA method.

    :return: A tuple with the index of the subfield of Fx, the starting
        position in it, the index of the subfield of Fy, the starting position
        in it, and the length of the match; or `None` if no match is found.
    """

    max_len_y = max(len(sf_y) for sf_y in f_y)
    for idx_x, sf_x in enumerate(f_x):
        # Binary search for the longest window of `sf_x` found in Fy
        low, high = 0, min(len(sf_x), max_len_y)
        while low < high:
            mid = (low + high + 1) // 2
            if _has_common_window(sf_x, f_y, mid):
                low = mid
            else:
                high = mid - 1

        if low:
            # Return the first window (in order) and its first occurrence
            for i in range(len(sf_x) - low + 1):
                pattern = sf_x[i : i + low]
                for idx_y, sf_y in enumerate(f_y):
                    j = sf_y.find(pattern)
                    if j >= 0:
                        return idx_x, i, idx_y, j, low

    return None


def _has_common_window(sf_x: str, f_y: List[str], length: int) -> bool:
    """
    Checks whether any window of `sf_x` of a given length is found in Fy.
    """

    windows = {
        sf_y[j : j + length] for sf_y in f_y for j in range(len(sf_y) - length + 1)
    }

    return any(sf_x[i : i + length] in windows for i in range(len(sf_x) - length + 1))
