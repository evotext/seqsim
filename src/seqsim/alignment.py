"""
Module implementing methods for sequence comparison based on alignments.

The methods in this module align two sequences of arbitrary hashable elements
with user-configurable costs or scores, as commonly done in bioinformatics:
global alignment with affine gaps (Needleman-Wunsch, with the algorithm by
Gotoh), local alignment (Smith-Waterman), and the Monge-Elkan similarity for
sequences of sequences. See the `edit` module for the naming convention of
the functions.
"""

# Import Python standard libraries
from typing import Callable, Hashable, List, Sequence

# Import local modules
from .common import empty_dissim

_INF = float("inf")


def _unit_cost(elem_x: Hashable, elem_y: Hashable) -> float:
    """
    Default substitution cost: zero for equal elements, one otherwise.
    """

    return 0.0 if elem_x == elem_y else 1.0


def _unit_score(elem_x: Hashable, elem_y: Hashable) -> float:
    """
    Default substitution score: one for equal elements, minus one otherwise.
    """

    return 1.0 if elem_x == elem_y else -1.0


def _check_gaps(gap_open: float, gap_extend: float) -> None:
    """
    Raises a `ValueError` if gap costs are invalid.
    """

    if gap_open < 0:
        raise ValueError(f"`gap_open` must be non-negative, got {gap_open!r}.")
    if gap_extend <= 0:
        raise ValueError(f"`gap_extend` must be positive, got {gap_extend!r}.")


def nw_dissim(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    *,
    sub_cost: Callable[[Hashable, Hashable], float] = _unit_cost,
    gap_open: float = 0.0,
    gap_extend: float = 1.0,
    normal: bool = False,
) -> float:
    """
    Computes the cost of the optimal global alignment of two sequences.

    This is the Needleman-Wunsch global alignment, expressed as a cost to be
    minimized, with affine gaps computed with the algorithm by Gotoh (1982):
    aligning element `a` with element `b` costs `sub_cost(a, b)`, and a gap
    of length `k` (a run of consecutive insertions or deletions) costs
    `gap_open + k * gap_extend`. With the default costs (unit substitutions,
    `gap_open=0`, `gap_extend=1`), the result is the Levenshtein distance.

    `sub_cost` must be symmetric and non-negative, and return zero for equal
    elements; under these conditions the result is symmetric and zero for
    identical sequences. It is a true distance when `sub_cost` is itself a
    metric and `gap_open` is zero, but not in general (for example, affine
    gaps violate the triangle inequality), hence the `_dissim` name.

    Example
    ********

    .. code-block:: python

        >>> seqsim.alignment.nw_dissim("kitten", "sitting")
        3.0
        >>> seqsim.alignment.nw_dissim("abcdef", "abef", gap_open=2.0)
        4.0

    References
    ***********

    Needleman, Saul B.; Wunsch, Christian D. (1970). "A general method applicable to
    the search for similarities in the amino acid sequence of two proteins". Journal
    of Molecular Biology 48 (3): 443–53. doi:10.1016/0022-2836(70)90057-4

    Gotoh, Osamu (1982). "An improved algorithm for matching biological sequences".
    Journal of Molecular Biology 162 (3): 705–708. doi:10.1016/0022-2836(82)90398-9

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param sub_cost: A function returning the cost of aligning two elements.
        Defaults to zero for equal elements and one otherwise.
    :param gap_open: The cost of opening a gap. Defaults to zero.
    :param gap_extend: The cost of each element in a gap. Must be positive.
        Defaults to one.
    :param normal: Whether to normalize the cost in range [0..1], by dividing
        it by the cost of aligning each sequence entirely against a gap (which
        is always an upper bound).
    :return: The cost of the optimal global alignment.
    """

    _check_gaps(gap_open, gap_extend)
    len_x, len_y = len(seq_x), len(seq_y)

    # `match[i][j]`: best cost ending with x_i aligned to y_j; `del_x[i][j]`:
    # best cost ending with x_i aligned to a gap; `ins_y[i][j]`: best cost
    # ending with y_j aligned to a gap
    match = [[_INF] * (len_y + 1) for _ in range(len_x + 1)]
    del_x = [[_INF] * (len_y + 1) for _ in range(len_x + 1)]
    ins_y = [[_INF] * (len_y + 1) for _ in range(len_x + 1)]
    match[0][0] = 0.0
    for i in range(1, len_x + 1):
        del_x[i][0] = gap_open + i * gap_extend
    for j in range(1, len_y + 1):
        ins_y[0][j] = gap_open + j * gap_extend

    for i in range(1, len_x + 1):
        for j in range(1, len_y + 1):
            cost = sub_cost(seq_x[i - 1], seq_y[j - 1])
            if cost < 0:
                raise ValueError("`sub_cost` must return non-negative values.")
            match[i][j] = (
                min(match[i - 1][j - 1], del_x[i - 1][j - 1], ins_y[i - 1][j - 1])
                + cost
            )
            del_x[i][j] = min(
                match[i - 1][j] + gap_open + gap_extend,
                del_x[i - 1][j] + gap_extend,
                ins_y[i - 1][j] + gap_open + gap_extend,
            )
            ins_y[i][j] = min(
                match[i][j - 1] + gap_open + gap_extend,
                ins_y[i][j - 1] + gap_extend,
                del_x[i][j - 1] + gap_open + gap_extend,
            )

    dist = min(match[len_x][len_y], del_x[len_x][len_y], ins_y[len_x][len_y])

    if normal:
        bound = sum(
            gap_open + length * gap_extend for length in (len_x, len_y) if length
        )
        return dist / bound if bound else 0.0

    return float(dist)


def sw_simil(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    *,
    score: Callable[[Hashable, Hashable], float] = _unit_score,
    gap_open: float = 0.0,
    gap_extend: float = 1.0,
    normal: bool = False,
) -> float:
    """
    Computes the score of the optimal local alignment of two sequences.

    This is the Smith-Waterman local alignment, with affine gaps (Gotoh,
    1982): it finds the pair of sub-sequences with the highest alignment
    score, where aligning element `a` with element `b` scores `score(a, b)`,
    and a gap of length `k` is penalized by `gap_open + k * gap_extend`. It
    is useful for detecting a shared run of elements within longer sequences.
    `score` must be symmetric, so that the result is symmetric.

    Example
    ********

    .. code-block:: python

        >>> seqsim.alignment.sw_simil("XXXabcdYYY", "ZZabcdZZ")
        4.0
        >>> seqsim.alignment.sw_simil("XXXabcdYYY", "ZZabcdZZ", normal=True)
        0.4

    References
    ***********

    Smith, Temple F.; Waterman, Michael S. (1981). "Identification of Common
    Molecular Subsequences". Journal of Molecular Biology 147: 195–197.
    doi:10.1016/0022-2836(81)90087-5

    Gotoh, Osamu (1982). "An improved algorithm for matching biological sequences".
    Journal of Molecular Biology 162 (3): 705–708. doi:10.1016/0022-2836(82)90398-9

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param score: A function returning the score of aligning two elements.
        Defaults to one for equal elements and minus one otherwise.
    :param gap_open: The penalty for opening a gap. Defaults to zero.
    :param gap_extend: The penalty for each element in a gap. Must be
        positive. Defaults to one.
    :param normal: Whether to normalize the score in range [0..1], by dividing
        it by the highest of the scores of each sequence aligned with itself.
        A normalized score of 1.0 indicates identical sequences (with the
        default scores).
    :return: The score of the optimal local alignment.
    """

    _check_gaps(gap_open, gap_extend)

    simil = _sw_score(seq_x, seq_y, score, gap_open, gap_extend)

    if normal:
        if not seq_x and not seq_y:
            return 1.0
        bound = max(
            _sw_score(seq_x, seq_x, score, gap_open, gap_extend),
            _sw_score(seq_y, seq_y, score, gap_open, gap_extend),
        )
        return min(simil / bound, 1.0) if bound > 0 else 0.0

    return float(simil)


def _sw_score(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    score: Callable[[Hashable, Hashable], float],
    gap_open: float,
    gap_extend: float,
) -> float:
    """
    Computes the Smith-Waterman score with affine gaps.
    """

    len_y = len(seq_y)
    best = 0.0

    # Rolling rows: `h` holds the best local score ending at (i, j), `e` the
    # best score ending with a gap in `seq_x`, and `f` with a gap in `seq_y`
    prev_h: List[float] = [0.0] * (len_y + 1)
    prev_f: List[float] = [-_INF] * (len_y + 1)
    for elem_x in seq_x:
        curr_h = [0.0] * (len_y + 1)
        curr_f = [-_INF] * (len_y + 1)
        e = -_INF
        for j in range(1, len_y + 1):
            e = max(curr_h[j - 1] - gap_open - gap_extend, e - gap_extend)
            curr_f[j] = max(prev_h[j] - gap_open - gap_extend, prev_f[j] - gap_extend)
            curr_h[j] = max(
                0.0,
                prev_h[j - 1] + score(elem_x, seq_y[j - 1]),
                e,
                curr_f[j],
            )
            best = max(best, curr_h[j])
        prev_h, prev_f = curr_h, curr_f

    return best


def _levenshtein_simil(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
    """
    Default inner similarity for Monge-Elkan: one minus the normalized Levenshtein.
    """

    # Imported here to avoid a circular import at module level
    from .edit import levenshtein_dist

    return 1.0 - levenshtein_dist(seq_x, seq_y, normal=True)


def monge_elkan_simil(
    seq_x: Sequence[Sequence[Hashable]],
    seq_y: Sequence[Sequence[Hashable]],
    *,
    inner: Callable[[Sequence[Hashable], Sequence[Hashable]], float] = (
        _levenshtein_simil
    ),
    normal: bool = False,
) -> float:
    """
    Computes the Monge-Elkan similarity between two sequences of sequences.

    Each element of the sequences is itself a sequence (for example, the
    titles or incipits of the texts in a manuscript, each a string). The
    Monge-Elkan similarity of `x` to `y` is the mean, over the elements of
    `x`, of their highest `inner` similarity to any element of `y`. As this
    is not symmetric, the mean of both directions is returned. The order of
    the elements is not taken into account.

    `inner` must return similarities in range [0..1]; by default it is one
    minus the normalized Levenshtein distance. Results are always in range
    [0..1], so `normal` has no effect.

    Example
    ********

    .. code-block:: python

        >>> seqsim.alignment.monge_elkan_simil(["vita antonii"], ["vita antonij"])
        0.9166666666666666

    References
    ***********

    Monge, Alvaro E.; Elkan, Charles (1996). "The Field Matching Problem: Algorithms
    and Applications". Proceedings of the Second International Conference on
    Knowledge Discovery and Data Mining (KDD'96): 267–270.

    :param seq_x: The first sequence of sequences to be compared.
    :param seq_y: The second sequence of sequences to be compared.
    :param inner: The similarity function, in range [0..1], used to compare
        the elements.
    :param normal: Ignored, as results are always in range [0..1].
    :return: The symmetric Monge-Elkan similarity.
    """

    empty = empty_dissim(seq_x, seq_y)
    if empty is not None:
        return 1.0 - empty

    def directional(source, target):
        return sum(max(inner(a, b) for b in target) for a in source) / len(source)

    return (directional(seq_x, seq_y) + directional(seq_y, seq_x)) / 2.0
