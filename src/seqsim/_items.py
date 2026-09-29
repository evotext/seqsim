"""
Views of a sequence as items: occurrences, shared items, adjacencies, n-grams.

The measures of order and of tokens, and the analysis of traditions, do not
compare sequences element by element, but through views of their items. This
module computes these views in one place (see `CONTEXT.md` for the
vocabulary):

  * occurrence labels, distinguishing repeated elements (`("a", 0)` for the
    first `"a"`, `("a", 1)` for the second, and so on);
  * first occurrences, reducing a sequence to its distinct items;
  * the restriction of two sequences to the items they share;
  * adjacencies, the pairs of consecutive elements, optionally with the
    boundaries of the sequence;
  * n-grams, optionally padded with the boundary.
"""

# Import Python standard libraries
from collections import Counter
from typing import Hashable, Iterable, List, Sequence, Tuple

Label = Tuple[Hashable, int]

REPEATS = ("occurrence", "first")


class _Boundary:
    """
    Sentinel marking the start and the end of a sequence.

    A dedicated object is used, instead of a string such as `"$$$"`, so that
    a boundary can never be confused with an element of a sequence. There is
    a single instance, which survives pickling, exported as `ngrams.PAD`.
    """

    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __repr__(self) -> str:
        return "PAD"

    def __reduce__(self):
        return (_Boundary, ())


BOUNDARY = _Boundary()


def occurrences(seq: Iterable[Hashable]) -> List[Label]:
    """
    Labels each element with its occurrence number (0 for the first).
    """

    counter: Counter = Counter()
    labels = []
    for element in seq:
        labels.append((element, counter[element]))
        counter[element] += 1

    return labels


def first_occurrences(seq: Iterable[Hashable]) -> List[Hashable]:
    """
    Returns the distinct elements, in the order of their first occurrence.
    """

    return list(dict.fromkeys(seq))


def shared(
    seq_x: Sequence[Hashable],
    seq_y: Sequence[Hashable],
    *,
    repeats: str = "occurrence",
) -> Tuple[List[Label], List[Label]]:
    """
    Restricts two sequences to the occurrences they share, each in its order.

    With `repeats="first"`, both sequences are first reduced to the first
    occurrence of each element.

    :return: The occurrence labels of the shared elements of each sequence.
    """

    if repeats == "first":
        seq_x, seq_y = first_occurrences(seq_x), first_occurrences(seq_y)
    elif repeats not in REPEATS:
        raise ValueError(f"`repeats` must be 'occurrence' or 'first', got {repeats!r}.")

    labels_x, labels_y = occurrences(seq_x), occurrences(seq_y)
    common = set(labels_x) & set(labels_y)

    return (
        [label for label in labels_x if label in common],
        [label for label in labels_y if label in common],
    )


def adjacencies(seq: Iterable[Hashable], *, boundaries: bool) -> Counter:
    """
    Returns the multiset of ordered pairs of consecutive elements.

    With `boundaries`, the sequence is framed by `BOUNDARY`, so that a
    sequence of `n` elements has `n + 1` adjacencies.
    """

    items = [BOUNDARY, *seq, BOUNDARY] if boundaries else list(seq)

    return Counter(zip(items, items[1:]))


def ngrams(seq: Iterable[Hashable], size: int, *, pad: bool = False) -> Counter:
    """
    Returns the multiset of contiguous sub-sequences of `size` elements.

    With `pad`, the sequence is framed by `size - 1` boundaries at each end,
    so that every element appears in `size` n-grams.
    """

    items = tuple(seq)
    if pad:
        items = (BOUNDARY,) * (size - 1) + items + (BOUNDARY,) * (size - 1)

    return Counter(items[idx : idx + size] for idx in range(len(items) - size + 1))
