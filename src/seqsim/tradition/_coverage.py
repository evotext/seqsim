"""
Coverage: where the absence of an item from a witness is evidence.

An item absent from a witness is evidence about transmission only if the
witness covers the place where the item belongs. An absence inside a lost
quire, in a section the witness never had, or before a mutilated beginning is
missing data (`None`), not an absence (`0`).

Coverage is defined in a reference frame (see `consensus_order` and `reference_labels`): a
consensus order of all items, and a label (for example a chapter) for each of
them. A witness covers a label if it has at least `min_share` of the items
with that label, and at least `min_items` of them. Outside covered labels,
before the first and after the last item of the witness in frame order, and
inside a recorded gap (a lacuna), an absence is missing data. The share
threshold also turns absences into missing data for excerpt collections,
which take a few items from many sections: they cover no section densely
enough for an absence to mean anything.
"""

# Import Python standard libraries
from collections import Counter
from typing import (
    Callable,
    Dict,
    Hashable,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
    Tuple,
)

# Import local modules
from ._frame import unique

# A cell of a character matrix: 1 (present), 0 (absent), or None (missing)
Cell = Optional[int]


def covered_labels(
    sequence: Sequence[Hashable],
    frame: Sequence[Hashable],
    label_of: Mapping[Hashable, Hashable],
    min_items: int = 2,
    min_share: float = 0.06,
) -> Set[Hashable]:
    """
    Returns the labels a witness covers.

    A witness covers a label if it has at least `min_share` of the items of
    the frame with that label, and at least `min_items` of them. Items
    without a label (and the label `None`) are never covered.

    Example
    ********

    .. code-block:: python

        >>> label_of = {"a": "I", "b": "I", "c": "II", "d": "II", "e": "II"}
        >>> sorted(seqsim.tradition.covered_labels("abc", "abcde", label_of))
        ['I']
        >>> sorted(seqsim.tradition.covered_labels("abc", "abcde", label_of, min_items=1))
        ['I', 'II']

    :param sequence: The items of the witness.
    :param frame: The items of the reference frame.
    :param label_of: A mapping from items to their label.
    :param min_items: The minimum number of items of a label. Defaults to 2.
    :param min_share: The minimum share of the items of a label. Defaults to
        0.06.
    :return: The set of covered labels.
    """

    label_size = Counter(label_of.get(item) for item in frame)
    in_frame = set(frame)
    per_label = Counter(
        label_of.get(item) for item in set(sequence) if item in in_frame
    )

    return {
        label
        for label, count in per_label.items()
        if label is not None
        and count >= min_items
        and count >= min_share * label_size[label]
    }


def coverage(
    sequence: Sequence[Hashable],
    frame: Sequence[Hashable],
    label_of: Mapping[Hashable, Hashable],
    gaps: Iterable[int] = (),
    min_items: int = 2,
    min_share: float = 0.06,
) -> Dict[Hashable, Cell]:
    """
    Codes every item of the frame as present, absent, or missing in a witness.

    Every item of `frame` gets a value: 1 if the witness has it; otherwise
    `None` (missing data) if its label is one the witness does not cover
    (see `covered_labels`), if it lies before the first or after the last
    item of the witness in frame order, or if it lies in the frame-order gap
    where a lacuna of the witness falls; and 0 (absent) otherwise.

    The witness is first reduced to the first occurrence of each item. A
    lacuna is given by the index, in this reduced sequence, of the item
    before which it falls (0 for the start, its length for the end), as
    computed by `proportional_positions` with the length of the reduced
    sequence. The frame items between the last item of the witness before
    the lacuna and the first one after it, in frame order, are missing.

    Example
    ********

    .. code-block:: python

        >>> frame = "abcdefgh"
        >>> label_of = dict(zip(frame, "IIIIIIII"))
        >>> seqsim.tradition.coverage("bdeg", frame, label_of)
        {'a': None, 'b': 1, 'c': 0, 'd': 1, 'e': 1, 'f': 0, 'g': 1, 'h': None}
        >>> seqsim.tradition.coverage("bdeg", frame, label_of, gaps=[3])
        {'a': None, 'b': 1, 'c': 0, 'd': 1, 'e': 1, 'f': None, 'g': 1, 'h': None}

    :param sequence: The items of the witness, in its order.
    :param frame: The items of the reference frame, in frame order.
    :param label_of: A mapping from items to their label.
    :param gaps: The indices, in the witness reduced to first occurrences, of
        the items before which a lacuna falls.
    :param min_items: The minimum number of items of a covered label.
        Defaults to 2.
    :param min_share: The minimum share of the items of a covered label.
        Defaults to 0.06.
    :return: A dictionary from each item of the frame, in frame order, to 1,
        0, or `None`.
    """

    order = unique(sequence)
    position = {item: pos for pos, item in enumerate(frame)}
    present = [item for item in order if item in position]
    if not present:
        return dict.fromkeys(frame)

    covered = covered_labels(order, frame, label_of, min_items, min_share)
    first = min(position[item] for item in present)
    last = max(position[item] for item in present)

    lost: Set[int] = set()
    for cut in gaps:
        before = [position[item] for item in order[:cut] if item in position]
        after = [position[item] for item in order[cut:] if item in position]
        low = before[-1] if before else -1
        high = after[0] if after else len(frame)
        if low < high:
            lost.update(range(low + 1, high))

    content = set(order)
    result: Dict[Hashable, Cell] = {}
    for pos, item in enumerate(frame):
        if item in content:
            result[item] = 1
        elif (
            label_of.get(item) not in covered
            or pos < first
            or pos > last
            or pos in lost
        ):
            result[item] = None
        else:
            result[item] = 0

    return result


def merge_coverage(
    coverages: Sequence[Mapping[Hashable, Cell]], frame: Sequence[Hashable]
) -> Dict[Hashable, Cell]:
    """
    Combines the coverage of the parts of a single witness.

    For a manuscript made of several parts (for example collections bound
    together, or sections catalogued separately), an item is present if any
    part has it, absent if no part has it and some part covers its place,
    and missing otherwise.

    Example
    ********

    .. code-block:: python

        >>> part_1 = {"a": 1, "b": 0, "c": None}
        >>> part_2 = {"a": None, "b": None, "c": 1}
        >>> seqsim.tradition.merge_coverage([part_1, part_2], "abc")
        {'a': 1, 'b': 0, 'c': 1}

    :param coverages: The coverage of each part, as returned by `coverage`.
    :param frame: The items of the reference frame.
    :return: A dictionary from each item of the frame to 1, 0, or `None`.
        With a single part, a copy of its coverage is returned.
    """

    if len(coverages) == 1:
        return dict(coverages[0])

    combined: Dict[Hashable, Cell] = {}
    for item in frame:
        values = [cells[item] for cells in coverages]
        combined[item] = 1 if 1 in values else 0 if 0 in values else None

    return combined


def density(coverage: Mapping[Hashable, Cell]) -> float:
    """
    Returns the share of the items a witness has among those it covers.

    The number of items present divided by the number of items present or
    absent, ignoring missing data; 0.0 if there are none.

    Example
    ********

    .. code-block:: python

        >>> seqsim.tradition.density({"a": 1, "b": 0, "c": None, "d": 1})
        0.6666666666666666

    :param coverage: A coverage, as returned by `coverage`.
    :return: The density of the witness.
    """

    values = list(coverage.values())
    present, absent = values.count(1), values.count(0)

    return present / (present + absent) if present + absent else 0.0


def merge_parts(
    parts: Sequence[Tuple[str, Sequence[Hashable]]],
    max_overlap: float = 0.2,
    compatible: Optional[Callable[[str, str], bool]] = None,
) -> Tuple[List[Tuple[str, Sequence[Hashable]]], List[Tuple[str, Sequence[Hashable]]]]:
    """
    Merges the complementary parts of a single codex into one witness.

    Catalogues often split one collection copied in a codex into
    consecutive or interleaved parts. Parts are considered from the largest
    (by number of distinct items, so that repetitions do not count; ties keep
    the input order): the largest starts the
    merged witness, and each following part is merged if less than
    `max_overlap` of its distinct items are already in the merged content,
    and if `compatible(part_label, largest_label)` holds (for example, if
    both are catalogued as the same kind of collection). Other parts, mostly
    contained in the merged content, are second copies of a section,
    possibly from another exemplar, and stay separate witnesses. An empty
    part has no overlap with the merged content, so it is merged (if
    `compatible` allows it) and its label appears in the merged label; this
    is intentional, and avoids the division by zero of a naive overlap ratio.

    The merged witness is made of its parts in part-label order: its label
    is usually their labels joined (e.g. `"A+C"`), its items their sequences
    concatenated in that order, and its coverage the combination of theirs
    (see `merge_coverage`).

    Example
    ********

    .. code-block:: python

        >>> parts = [("A", "abcd"), ("B", "efg"), ("C", "abx")]
        >>> merged, separate = seqsim.tradition.merge_parts(parts)
        >>> merged
        [('A', 'abcd'), ('B', 'efg')]
        >>> separate
        [('C', 'abx')]

    :param parts: The parts of the codex, as `(label, items)` pairs, with
        string labels.
    :param max_overlap: The maximum share of the distinct items of a part
        already in the merged content for merging it. Defaults to 0.2.
    :param compatible: An optional function of two part labels, telling
        whether a part can be merged with the largest part.
    :return: A tuple with the parts merged, in label order (an empty list if
        no part was merged with the largest one), and the parts left
        separate (all parts, in the input order, if none was merged; else in
        order of decreasing size).
    """

    parts = list(parts)
    if len(parts) < 2:
        return [], parts

    # Parts are sized by their number of distinct items, so that a part with
    # many repeated items does not become the largest one; the sort is
    # stable, so ties keep the input order
    by_size = sorted(parts, key=lambda part: len(set(part[1])), reverse=True)
    merged = [by_size[0]]
    content = set(by_size[0][1])
    separate = []
    for label, items in by_size[1:]:
        items_set = set(items)
        overlap = len(items_set & content) / len(items_set) if items_set else 0.0
        if overlap < max_overlap and (
            compatible is None or compatible(label, merged[0][0])
        ):
            merged.append((label, items))
            content |= items_set
        else:
            separate.append((label, items))

    if len(merged) == 1:
        return [], parts

    merged.sort(key=lambda part: part[0] or "")

    return merged, separate
