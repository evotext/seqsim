"""
Module implementing methods for sequence dissimilarity based on tokens.

These methods compare the elements (or sub-sequences) that two sequences
share, such as the Jaccard index, operating on arbitrary sequences of hashable
elements. See the `edit` module for the naming convention of the functions.
"""

# Import Python standard libraries
from collections import Counter
from typing import Hashable, Sequence

# Import local modules
from ._measure import Scored, measure
from .common import equivalent_string
from .ngrams import PAD


def _check_shingle_size(size: int) -> None:
    """
    Raises a `ValueError` if a shingle (q-gram) size is not a positive integer.
    """

    if isinstance(size, bool) or not isinstance(size, int) or size < 1:
        raise ValueError(f"The size must be a positive integer, got {size!r}.")


def _shingles(seq: Sequence[Hashable], size: int) -> Counter:
    """
    Returns the multiset of contiguous sub-sequences of a given size.
    """

    items = tuple(seq)

    return Counter(items[i : i + size] for i in range(len(items) - size + 1))


def _check_qgram_options(q: int, pad: bool) -> None:
    """
    Validates the options of `qgram_dissim`.
    """

    _check_shingle_size(q)


def _check_containment_options(size: int) -> None:
    """
    Validates the options of `containment`.
    """

    _check_shingle_size(size)


@measure(
    key="jaccard",
    kind="dissim",
    identity="no",
    identity_example=("a", "aa"),
    triangle="yes",
    empty="max",
)
def jaccard_dissim(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
    """
    Computes the Jaccard dissimilarity between two sequences.

    The dissimilarity is one minus the Jaccard index of the sets of elements
    of both sequences. While the Jaccard distance is a true distance on sets,
    order and repetition are ignored, so different sequences can have a
    dissimilarity of zero (e.g., `"ab"` and `"ba"`, or `"a"` and `"aa"`).

    Example
    ********

    .. code-block:: python

        >>> seqsim.token.jaccard_dissim("abc", "bcde")
        0.6

    References
    ***********

    Tan PN, Steinbach M, Kumar V (2005). Introduction to Data Mining. ISBN 0-321-32136-7.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :return: The Jaccard dissimilarity between the two sequences.
    """

    set_x, set_y = set(seq_x), set(seq_y)

    return 1.0 - (len(set_x & set_y) / len(set_x | set_y))


@measure(key="subseq_jaccard", kind="dissim", triangle="unproven", empty="max")
def subseq_jaccard_dissim(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]
) -> float:
    """
    Computes a Jaccard dissimilarity between two sequences using sub-sequences.

    For each length `n` from 1 to the length of the longest sequence, the
    Jaccard index is computed on the multisets of contiguous sub-sequences of
    length `n` of both sequences. The similarity is the mean of these indices,
    weighted by `n` so that longer shared sub-sequences count more, and the
    dissimilarity is one minus this similarity. Identical sequences, and only
    identical sequences, have a dissimilarity of zero; the measure does not
    satisfy the triangle inequality.

    Example
    ********

    .. code-block:: python

        >>> seqsim.token.subseq_jaccard_dissim("abc", "bcde")
        0.91

    References
    ***********

    Tan PN, Steinbach M, Kumar V (2005). Introduction to Data Mining. ISBN 0-321-32136-7.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :return: The Subseq-Jaccard dissimilarity between the two sequences.
    """

    # Strings are much faster to slice and hash than tuples
    str_x, str_y = equivalent_string(seq_x, seq_y)
    max_length = max(len(str_x), len(str_y))

    weighted_sum = 0.0
    for length in range(1, max_length + 1):
        counter_x = Counter(
            str_x[i : i + length] for i in range(len(str_x) - length + 1)
        )
        counter_y = Counter(
            str_y[i : i + length] for i in range(len(str_y) - length + 1)
        )

        # Use multisets for both the intersection and the union, so that
        # repeated sub-sequences are counted consistently
        intersection = sum((counter_x & counter_y).values())
        union = sum((counter_x | counter_y).values())
        weighted_sum += length * (intersection / union)

    # The highest possible value is the sum of all weights
    return 1.0 - (weighted_sum / ((max_length * (max_length + 1)) / 2))


@measure(
    key="sorensen",
    kind="dissim",
    identity="no",
    identity_example=("ab", "ba"),
    triangle="no",
    triangle_example=("cc", "aac", "aba"),
    empty="max",
)
def sorensen_dissim(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> float:
    """
    Computes a dissimilarity between two sequences based on the Sørensen–Dice coefficient.

    The dissimilarity is one minus the Sørensen–Dice coefficient of the
    multisets of elements of both sequences. Order is ignored, so different
    sequences can have a dissimilarity of zero (e.g., `"ab"` and `"ba"`), and
    the measure does not satisfy the triangle inequality.

    Example
    ********

    .. code-block:: python

        >>> seqsim.token.sorensen_dissim("abc", "bcde")
        0.4285714285714286

    References
    ***********

    Kondrak, Grzegorz; Marcu, Daniel; Knight, Kevin (2003). "Cognates Can Improve Statistical
    Translation Models" (PDF). Proceedings of HLT-NAACL 2003: Human Language Technology
    Conference of the North American Chapter of the Association for Computational
    Linguistics. pp. 46–48.

    Sørensen, T. (1948). "A method of establishing groups of equal amplitude in plant
    sociology based on similarity of species and its application to analyses of the
    vegetation on Danish commons". Kongelige Danske Videnskabernes Selskab. 5 (4): 1–34.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :return: The Sørensen–Dice dissimilarity between the two sequences.
    """

    intersection = sum((Counter(seq_x) & Counter(seq_y)).values())

    return 1.0 - (2.0 * intersection / (len(seq_x) + len(seq_y)))


@measure(
    key="qgram",
    kind="dissim",
    identity="no",
    identity_example=("abaca", "acaba"),
    triangle="yes",
    bound="scored",
    check=_check_qgram_options,
    normal_doc=(
        "Whether to normalize the dissimilarity in range [0..1] by dividing it "
        "by the total number of q-grams in both sequences."
    ),
)
def qgram_dissim(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    *,
    q: int = 2,
    pad: bool = True,
) -> float:
    """
    Computes the q-gram dissimilarity between two sequences.

    The q-gram distance of Ukkonen (1992) is the L1 distance between the
    q-gram profiles of both sequences, that is, the sum over all contiguous
    sub-sequences of `q` elements of the absolute difference of their number
    of occurrences in each sequence. It is a lower bound for the edit
    distance and can be computed in linear time. By default, sequences are
    padded with `q - 1` boundary symbols on each side (as in the `ngrams`
    module), so that elements at the boundaries are counted as often as the
    others and every non-empty sequence has q-grams; `pad=False` gives
    Ukkonen's original definition. Different sequences can have the same
    profile (e.g., `"abaca"` and `"acaba"`), so identity of
    indiscernibles does not hold; the triangle inequality does.

    Example
    ********

    .. code-block:: python

        >>> seqsim.token.qgram_dissim("abc", "bca")
        6.0
        >>> seqsim.token.qgram_dissim("abc", "bca", pad=False)
        2.0

    References
    ***********

    Ukkonen, Esko (1992). "Approximate string-matching with q-grams and maximal
    matches". Theoretical Computer Science 92 (1): 191–211.
    doi:10.1016/0304-3975(92)90143-4

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param q: The number of elements in each q-gram. Defaults to 2.
    :param pad: Whether to pad the sequences with boundary symbols. Defaults
        to `True`.
    :return: The q-gram dissimilarity.
    """

    if pad:
        seq_x = [PAD] * (q - 1) + list(seq_x) + [PAD] * (q - 1)
        seq_y = [PAD] * (q - 1) + list(seq_y) + [PAD] * (q - 1)
    grams_x, grams_y = _shingles(seq_x, q), _shingles(seq_y, q)

    dist = sum((grams_x - grams_y).values()) + sum((grams_y - grams_x).values())
    total = sum(grams_x.values()) + sum(grams_y.values())

    return Scored(dist, total)


@measure(kind="simil", symmetric=False)
def tversky_simil(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    *,
    alpha: float = 0.5,
    beta: float = 0.5,
) -> float:
    """
    Computes the Tversky index between two sequences.

    The Tversky index generalizes the Jaccard and Sørensen–Dice coefficients,
    weighting the elements unique to each sequence differently. On the
    multisets `X` and `Y` of elements it is `|X & Y| / (|X & Y| + alpha *
    |X - Y| + beta * |Y - X|)`. With `alpha = beta = 0.5` (the default) it is
    the Sørensen–Dice coefficient, and with `alpha = beta = 1` the (multiset)
    Jaccard index. It is symmetric only when `alpha == beta`: for example,
    with `alpha=1` and `beta=0` it measures how much of `x` is found in `y`.

    Example
    ********

    .. code-block:: python

        >>> seqsim.token.tversky_simil("abc", "abcdef", alpha=1.0, beta=0.0)
        1.0
        >>> seqsim.token.tversky_simil("abc", "abcdef")
        0.6666666666666666

    References
    ***********

    Tversky, Amos (1977). "Features of similarity". Psychological Review 84 (4):
    327–352. doi:10.1037/0033-295X.84.4.327

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param alpha: The weight of the elements unique to `seq_x`. Must be
        non-negative.
    :param beta: The weight of the elements unique to `seq_y`. Must be
        non-negative.
    :return: The Tversky index.
    """

    if alpha < 0 or beta < 0:
        raise ValueError("`alpha` and `beta` must be non-negative.")

    counter_x, counter_y = Counter(seq_x), Counter(seq_y)
    common = sum((counter_x & counter_y).values())
    only_x = sum((counter_x - counter_y).values())
    only_y = sum((counter_y - counter_x).values())

    denominator = common + alpha * only_x + beta * only_y
    if denominator == 0:
        # Both empty, or no common element with zero weights
        return 1.0 if not only_x and not only_y else 0.0

    return common / denominator


# Identical sequences have a containment of 1.0, the largest value
@measure(
    kind="directional",
    identical_zero=False,
    symmetric=False,
    normal=False,
    check=_check_containment_options,
)
def containment(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    *,
    size: int = 1,
) -> float:
    """
    Computes how much of `seq_x` is contained in `seq_y`.

    This is the containment of Broder (1997): the proportion of the
    contiguous sub-sequences ("shingles") of `size` elements of `seq_x` that
    are also found in `seq_y`, counted as multisets. It is directional by
    design, answering questions such as "is manuscript `x` an excerpt of
    manuscript `y`?", and thus not part of the naming convention for
    dissimilarities. An empty sequence (or one shorter than `size`) is
    contained in any sequence.

    Example
    ********

    .. code-block:: python

        >>> seqsim.token.containment("bcd", "abcdef", size=2)
        1.0
        >>> seqsim.token.containment("abcdef", "bcd", size=2)
        0.4

    References
    ***********

    Broder, Andrei Z. (1997). "On the resemblance and containment of documents".
    Proceedings of Compression and Complexity of SEQUENCES 1997: 21–29.
    doi:10.1109/SEQUEN.1997.666900

    :param seq_x: The sequence whose containment is measured.
    :param seq_y: The sequence in which `seq_x` is searched.
    :param size: The number of elements in each shingle. Defaults to 1.
    :return: The containment of `seq_x` in `seq_y`, in range [0..1].
    """

    shingles_x, shingles_y = _shingles(seq_x, size), _shingles(seq_y, size)
    total = sum(shingles_x.values())
    if not total:
        return 1.0

    return sum((shingles_x & shingles_y).values()) / total
