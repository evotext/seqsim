"""
Module implementing methods for comparing the order of elements in sequences.

These methods come from the comparison of rankings (e.g., Kendall's tau) and
of genome rearrangements (e.g., breakpoints), and are suited to comparing the
order of shared items, such as the texts in manuscripts with overlapping
contents. They are classically defined for permutations of the same distinct
items; here they are generalized to arbitrary sequences as follows:

  * repeated elements are distinguished by their occurrence, so that the
    first `"a"` of one sequence corresponds to the first `"a"` of the other,
    the second to the second, and so on;
  * elements (or occurrences) present in only one of the sequences are
    handled as described in each function.

Matching occurrences in order is a heuristic: for some measures, the optimal
correspondence between repeated elements is computationally hard to find.

See the `edit` module for the naming convention of the functions.
"""

# Import Python standard libraries
from bisect import bisect_left
from collections import Counter
from typing import Dict, Hashable, List, Optional, Sequence, Tuple

# Import local modules
from .common import lcs_length


class _Boundary:
    """
    Sentinel for sequence boundaries in adjacencies.
    """

    def __init__(self, name: str):
        self.name = name

    def __repr__(self) -> str:
        return self.name


_START = _Boundary("START")
_END = _Boundary("END")


def _occurrences(seq: Sequence[Hashable]) -> List[Tuple[Hashable, int]]:
    """
    Labels each element with its occurrence number (0 for the first).
    """

    counter: Counter = Counter()
    labels = []
    for element in seq:
        labels.append((element, counter[element]))
        counter[element] += 1

    return labels


def _shared_permutation(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]
) -> Tuple[List[int], int]:
    """
    Maps the shared (labelled) elements of `seq_x` to their positions in `seq_y`.

    :return: A tuple with the list of positions in `seq_y` (restricted to
        shared elements, renumbered from zero) of the shared elements in the
        order of `seq_x`, and the number of unshared elements in both
        sequences.
    """

    labels_x, labels_y = _occurrences(seq_x), _occurrences(seq_y)
    shared = set(labels_x) & set(labels_y)
    position_y: Dict[Tuple[Hashable, int], int] = {
        label: pos
        for pos, label in enumerate(label for label in labels_y if label in shared)
    }
    perm = [position_y[label] for label in labels_x if label in shared]
    unshared = len(labels_x) + len(labels_y) - 2 * len(shared)

    return perm, unshared


def _normalize(dist: float, bound: float, normal: bool) -> float:
    """
    Returns the distance as a float, divided by `bound` if requested.
    """

    if normal:
        return dist / bound if bound else 0.0

    return float(dist)


def kendall_tau_dissim(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    *,
    p: float = 0.5,
    normal: bool = False,
) -> float:
    """
    Computes the Kendall tau distance between two sequences.

    For two permutations of the same items, the Kendall tau distance is the
    number of pairs of items in a different relative order, which is also
    the minimum number of swaps of adjacent items needed to transform one
    into the other. For sequences with different items, it follows the
    generalization of Fagin et al. (2003) for top-k lists, where items
    missing from a sequence are considered to be ranked after all of its
    items. For each pair of items:

      * both items in both sequences: 1 if their order differs;
      * both items in one sequence, only one of them in the other: 1 if, in
        the first sequence, the item missing from the other comes first;
      * one item only in each sequence: 1;
      * both items in one sequence and none in the other: `p`, as their
        relative order in the other sequence is unknown.

    An end boundary, present in both sequences, is included as an item, so
    that an item present in only one sequence always counts. For two
    permutations of the same items this is a true distance, but the
    generalization does not satisfy the triangle inequality (e.g., `"ab"`,
    `"ac"`, and `"cd"`), hence the `_dissim` name. `p=0.5` is the "neutral"
    choice of Fagin et al.

    Example
    ********

    .. code-block:: python

        >>> seqsim.order.kendall_tau_dissim("abcd", "bacd")
        1.0
        >>> seqsim.order.kendall_tau_dissim("abcd", "dcba")
        6.0

    References
    ***********

    Kendall, Maurice G. (1938). "A New Measure of Rank Correlation". Biometrika 30
    (1–2): 81–89.

    Fagin, Ronald; Kumar, Ravi; Sivakumar, D. (2003). "Comparing top k lists".
    SIAM Journal on Discrete Mathematics 17 (1): 134–160.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param p: The penalty for pairs of items whose relative order is unknown,
        in range [0..1]. Defaults to 0.5.
    :param normal: Whether to normalize the distance in range [0..1] by
        dividing it by the number of pairs of distinct (labelled) items.
    :return: The Kendall tau dissimilarity.
    """

    if not 0.0 <= p <= 1.0:
        raise ValueError(f"`p` must be in range [0..1], got {p!r}.")

    # An end boundary, present in both sequences, makes items present in only
    # one sequence count even when there are no other items
    rank_x = {label: pos for pos, label in enumerate([*_occurrences(seq_x), _END])}
    rank_y = {label: pos for pos, label in enumerate([*_occurrences(seq_y), _END])}
    items = list(dict.fromkeys([*rank_x, *rank_y]))

    dist = 0.0
    for idx, item_i in enumerate(items):
        for item_j in items[idx + 1 :]:
            dist += _kendall_penalty(item_i, item_j, rank_x, rank_y, p)

    pairs = len(items) * (len(items) - 1) / 2
    return _normalize(dist, pairs, normal)


def _kendall_penalty(item_i, item_j, rank_x, rank_y, p) -> float:
    """
    Returns the penalty of a pair of items for the generalized Kendall tau.
    """

    in_x = (item_i in rank_x, item_j in rank_x)
    in_y = (item_i in rank_y, item_j in rank_y)

    # Both items in both sequences
    if all(in_x) and all(in_y):
        order_x = rank_x[item_i] < rank_x[item_j]
        order_y = rank_y[item_i] < rank_y[item_j]
        return float(order_x != order_y)

    # Both items in one sequence, exactly one of them in the other
    for ranks, other in ((rank_x, rank_y), (rank_y, rank_x)):
        if (
            item_i in ranks
            and item_j in ranks
            and (item_i in other) != (item_j in other)
        ):
            present = item_i if item_i in other else item_j
            absent = item_j if present is item_i else item_i
            return float(ranks[absent] < ranks[present])

    # Both items in one sequence, none in the other
    if all(in_x) or all(in_y):
        return p

    # One item only in each sequence
    return 1.0


def footrule_dissim(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    *,
    ell: Optional[int] = None,
    normal: bool = False,
) -> float:
    """
    Computes the Spearman footrule distance between two sequences.

    For two permutations of the same items, the footrule distance is the sum,
    over all items, of the absolute difference of their positions in both
    sequences. For sequences with different items, it follows the
    generalization of Fagin et al. (2003) for top-k lists: an item missing
    from a sequence is placed at position `ell` (counting from 1), by default
    one after the end of the longest sequence. The footrule is within a
    factor of two of the Kendall tau distance (Diaconis & Graham, 1977).

    With a fixed `ell`, larger than the length of all the sequences being
    compared, this is a true distance (an L1 distance between position
    vectors). With the default, `ell` depends on the pair of sequences and
    the triangle inequality can fail (e.g., `""`, `"a"`, and `"aa"`), hence
    the `_dissim` name.

    Example
    ********

    .. code-block:: python

        >>> seqsim.order.footrule_dissim("abcd", "bacd")
        2.0
        >>> seqsim.order.footrule_dissim("abcd", "dcba")
        8.0

    References
    ***********

    Diaconis, Persi; Graham, Ronald L. (1977). "Spearman's footrule as a measure of
    disarray". Journal of the Royal Statistical Society, Series B 39 (2): 262–268.
    doi:10.1111/j.2517-6161.1977.tb01624.x

    Fagin, Ronald; Kumar, Ravi; Sivakumar, D. (2003). "Comparing top k lists".
    SIAM Journal on Discrete Mathematics 17 (1): 134–160.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param ell: The position assigned to missing items. Must be larger than
        the length of both sequences. Defaults to one plus the length of the
        longest sequence.
    :param normal: Whether to normalize the distance in range [0..1] by
        dividing it by its value for two sequences with no item in common,
        which is its maximum.
    :return: The Spearman footrule dissimilarity.
    """

    max_len = max(len(seq_x), len(seq_y))
    if ell is None:
        ell = max_len + 1
    elif ell <= max_len:
        raise ValueError(
            f"`ell` must be larger than the length of both sequences, got {ell!r}."
        )

    rank_x = {label: pos for pos, label in enumerate(_occurrences(seq_x), start=1)}
    rank_y = {label: pos for pos, label in enumerate(_occurrences(seq_y), start=1)}
    items = set(rank_x) | set(rank_y)

    dist = sum(abs(rank_x.get(item, ell) - rank_y.get(item, ell)) for item in items)

    bound = sum(ell - pos for pos in rank_x.values()) + sum(
        ell - pos for pos in rank_y.values()
    )
    return _normalize(dist, bound, normal)


def ulam_dist(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], *, normal: bool = False
) -> float:
    """
    Computes the Ulam distance between two sequences.

    For two permutations of the same items, the Ulam distance is the minimum
    number of items that must be moved (taken out and reinserted elsewhere)
    to transform one into the other, which is the number of items minus the
    length of their longest common subsequence; it models, for example, a
    scribe moving single texts. Here it is generalized to the minimum number
    of moves, insertions, and deletions of single elements, which is
    `len(x) + len(y) - M - LCS(x, y)`, where `M` is the number of elements
    the sequences have in common (as multisets) and `LCS` the length of their
    longest common subsequence. As an edit distance with symmetric, unit-cost
    operations, it is a true distance.

    For sequences without repeated elements it is computed in
    `O(n log n)` time; otherwise, in `O(len(x) * len(y))`.

    Example
    ********

    .. code-block:: python

        >>> seqsim.order.ulam_dist("abcdef", "bcdefa")
        1.0
        >>> seqsim.order.ulam_dist("abcdef", "bcdXefa")
        2.0

    References
    ***********

    Aldous, David; Diaconis, Persi (1999). "Longest increasing subsequences: from
    patience sorting to the Baik-Deift-Johansson theorem". Bulletin of the American
    Mathematical Society 36 (4): 413–432.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Whether to normalize the distance in range [0..1] by
        dividing it by `len(x) + len(y) - M`, its maximum.
    :return: The Ulam distance.
    """

    common = sum((Counter(seq_x) & Counter(seq_y)).values())

    if len(set(seq_x)) == len(seq_x) and len(set(seq_y)) == len(seq_y):
        # Without repetitions, the LCS is the longest increasing subsequence
        # of the positions of the shared elements
        perm, _ = _shared_permutation(seq_x, seq_y)
        tails: List[int] = []
        for value in perm:
            idx = bisect_left(tails, value)
            if idx == len(tails):
                tails.append(value)
            else:
                tails[idx] = value
        lcs = len(tails)
    else:
        lcs = lcs_length(seq_x, seq_y)

    bound = len(seq_x) + len(seq_y) - common
    return _normalize(bound - lcs, bound, normal)


def cayley_dissim(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], *, normal: bool = False
) -> float:
    """
    Computes the Cayley distance between two sequences.

    For two permutations of the same items, the Cayley distance is the
    minimum number of swaps of any two items needed to transform one into
    the other, which is the number of items minus the number of cycles of
    the permutation relating them. Here, each element present in only one of
    the sequences adds one (for its deletion or insertion), and the Cayley
    distance is computed on the shared elements. This generalization does
    not satisfy the triangle inequality when elements are repeated (e.g.,
    `"aab"`, `"caab"`, and `"baac"`), hence the `_dissim` name; for
    permutations of the same items it is the Cayley distance, a true metric.

    Example
    ********

    .. code-block:: python

        >>> seqsim.order.cayley_dissim("abcdef", "fbcdea")
        1.0

    References
    ***********

    Cayley, Arthur (1849). "Note on the theory of permutations". Philosophical
    Magazine 34: 527–529.

    Diaconis, Persi (1988). Group Representations in Probability and Statistics.
    Institute of Mathematical Statistics, Lecture Notes 11.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Whether to normalize the distance in range [0..1] by
        dividing it by the number of distinct (labelled) elements in both
        sequences, which is an upper bound.
    :return: The Cayley dissimilarity.
    """

    perm, unshared = _shared_permutation(seq_x, seq_y)

    # Count the cycles of the permutation
    seen = [False] * len(perm)
    cycles = 0
    for start in range(len(perm)):
        if not seen[start]:
            cycles += 1
            pos = start
            while not seen[pos]:
                seen[pos] = True
                pos = perm[pos]

    dist = unshared + len(perm) - cycles

    return _normalize(dist, unshared + len(perm), normal)


def block_interchange_dissim(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], *, normal: bool = False
) -> float:
    """
    Computes the block interchange dissimilarity between two sequences.

    For two permutations of the same items, the block interchange distance is
    the minimum number of exchanges of two blocks of consecutive items (not
    necessarily adjacent, nor of the same size) needed to transform one into
    the other, modelling for example two groups of texts (or quires) swapped
    in a manuscript; it includes moving a single block elsewhere as a special
    case. It is computed with the cycle graph formula of Christie (1996),
    `(n + 1 - c) / 2`, where `c` is the number of cycles. Here, each element
    present in only one of the sequences adds one (for its deletion or
    insertion), and the distance is computed on the shared elements. This
    generalization does not always satisfy the triangle inequality, hence
    the `_dissim` name; for permutations of the same items it is a true
    distance.

    Example
    ********

    .. code-block:: python

        >>> seqsim.order.block_interchange_dissim("abcdefgh", "efghabcd")
        1.0
        >>> seqsim.order.block_interchange_dissim("abcdefgh", "agcdefbh")
        1.0

    References
    ***********

    Christie, David A. (1996). "Sorting permutations by block-interchanges".
    Information Processing Letters 60 (4): 165–169.
    doi:10.1016/S0020-0190(96)00155-X

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Whether to normalize the dissimilarity in range [0..1] by
        dividing it by the number of distinct (labelled) elements in both
        sequences, which is an upper bound.
    :return: The block interchange dissimilarity.
    """

    perm, unshared = _shared_permutation(seq_x, seq_y)
    size = len(perm)

    # Frame the permutation (with values from 1 to n) with 0 and n + 1, and
    # count the cycles of the cycle graph: from each value v, follow the
    # value preceding it in the permutation and move to its successor
    framed = [0, *[value + 1 for value in perm], size + 1]
    predecessor = {framed[i]: framed[i - 1] for i in range(1, size + 2)}
    seen = set()
    cycles = 0
    for start in range(1, size + 2):
        if start not in seen:
            cycles += 1
            value = start
            while value not in seen:
                seen.add(value)
                value = predecessor[value] + 1

    dist = unshared + (size + 1 - cycles) // 2

    return _normalize(dist, unshared + size, normal)


def breakpoint_dissim(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], *, normal: bool = False
) -> float:
    """
    Computes the breakpoint dissimilarity between two sequences.

    An adjacency is a pair of consecutive elements (`x[i]`, `x[i + 1]`),
    including a start and an end boundary, so that a sequence of length `n`
    has `n + 1` adjacencies. The dissimilarity is half the size of the
    symmetric difference of the multisets of adjacencies of both sequences.
    For two permutations of the same items, this is the classic breakpoint
    distance of genome rearrangements: the number of adjacencies of one that
    are broken in the other. Adjacencies are ordered (`"ab"` is different
    from `"ba"`), as the direction of reading matters for texts.

    When elements are repeated, different sequences can have the same
    adjacencies (e.g., `"abacada"` and `"acabada"`), so identity of
    indiscernibles does not hold; the triangle inequality does.

    Example
    ********

    .. code-block:: python

        >>> seqsim.order.breakpoint_dissim("abcdef", "abcfed")
        4.0

    References
    ***********

    Sankoff, David; Blanchette, Mathieu (1998). "Multiple genome rearrangement and
    breakpoint phylogeny". Journal of Computational Biology 5 (3): 555–570.
    doi:10.1089/cmb.1998.5.555

    Spencer, Matthew; Bordalejo, Barbara; Wang, Li-San; Barbrook, Adrian C.; Mooney,
    Linne R.; Robinson, Peter; Warnow, Tandy; Howe, Christopher J. (2003).
    "Analyzing the order of items in manuscripts of The Canterbury Tales". Computers
    and the Humanities 37 (1): 97–109. doi:10.1023/A:1021818600001

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Whether to normalize the dissimilarity in range [0..1] by
        dividing it by the mean number of adjacencies of both sequences.
    :return: The breakpoint dissimilarity.
    """

    adj_x, adj_y = _adjacencies(seq_x), _adjacencies(seq_y)
    diff = sum((adj_x - adj_y).values()) + sum((adj_y - adj_x).values())
    dist = diff / 2.0

    total = (sum(adj_x.values()) + sum(adj_y.values())) / 2.0
    return _normalize(dist, total, normal)


def _adjacencies(seq: Sequence[Hashable]) -> Counter:
    """
    Returns the multiset of (ordered) adjacencies of a sequence, with boundaries.
    """

    padded = [_START, *seq, _END]

    return Counter(zip(padded, padded[1:]))
