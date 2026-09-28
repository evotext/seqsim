"""
Reference frame of a tradition and the labels (chapters, books) of its items.

A tradition is a set of witnesses, each an ordered sequence of hashable
items (texts, sayings, strophes), that share a reference frame: a single
order of all the items, and a label for each of them, such as the chapter or
the book where it belongs. The frame and the labels decide where the absence
of an item from a witness is evidence (see `coverage`).

Witnesses are often described in another unit than their items: a catalogue
may count paragraphs, or folios, and give the chapter boundaries and lacunae
in that unit. The functions `proportional_labels` and
`proportional_positions` place such information in the sequence of items by
proportion, and `monotone_labels` and `reference_labels` refine the placement
with the agreement between witnesses.
"""

# Import Python standard libraries
from collections import Counter, defaultdict
from typing import (
    Dict,
    Hashable,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Tuple,
)


def unique(seq: Iterable[Hashable]) -> List[Hashable]:
    """
    Returns the items of a sequence in the order of their first occurrence.

    Example
    ********

    .. code-block:: python

        >>> seqsim.tradition.unique("abacb")
        ['a', 'b', 'c']

    :param seq: A sequence of hashable items.
    :return: The distinct items, in order of first occurrence.
    """

    return list(dict.fromkeys(seq))


def consensus_order(sequences: Iterable[Sequence[Hashable]]) -> List[Hashable]:
    """
    Merges the orders of several witnesses into one sequence of all their items.

    Witnesses are merged in the order given, each reduced to the first
    occurrence of its items. An item not yet placed goes directly after the
    nearest preceding item of the same witness that is already placed, or
    at the start if there is none. The result depends on the order of the
    witnesses: the first one fixes the order of its items, and the following
    ones only insert theirs. When the witnesses largely agree on the order,
    this is a good approximation to a consensus ranking; list first the
    witnesses whose order is most representative.

    Example
    ********

    .. code-block:: python

        >>> seqsim.tradition.consensus_order(["abde", "bcd", "aef"])
        ['a', 'b', 'c', 'd', 'e', 'f']

    :param sequences: The witnesses, as sequences of hashable items, in the
        order in which they are merged.
    :return: A list with all items of the witnesses, each once.
    """

    frame: List[Hashable] = []
    placed = set()
    for sequence in sequences:
        anchor = -1
        insertions: Dict[int, List[Hashable]] = defaultdict(list)
        index = {item: pos for pos, item in enumerate(frame)}
        for item in unique(sequence):
            if item in placed:
                anchor = index[item]
            else:
                insertions[anchor].append(item)
                placed.add(item)
        merged = list(insertions.get(-1, []))
        for pos, item in enumerate(frame):
            merged.append(item)
            merged.extend(insertions.get(pos, []))
        frame = merged

    return frame


def proportional_labels(
    n_items: int,
    total_units: float,
    events: Sequence[Tuple[float, Hashable]],
) -> List[Optional[Hashable]]:
    """
    Labels each item position from label boundaries given in another unit.

    A witness of `n_items` items is described in another unit (for example
    paragraphs of a catalogue), of which it has `total_units`, with each
    label (for example a chapter) starting at a given position in that
    unit. Item `i` is placed at position `i * total_units / n_items`, and
    gets the label of the last event, in the order given, whose position is
    lower than or equal to it; items before the first event get `None`.

    Example
    ********

    .. code-block:: python

        >>> seqsim.tradition.proportional_labels(6, 3, [(0, "I"), (2, "II")])
        ['I', 'I', 'I', 'I', 'II', 'II']
        >>> seqsim.tradition.proportional_labels(3, 3, [(1, "I")])
        [None, 'I', 'I']

    :param n_items: The number of items of the witness.
    :param total_units: The size of the witness in the other unit.
    :param events: A sequence of `(position, label)` pairs, in the other
        unit, giving where each label starts.
    :return: A list with the label of each item position, or `None`. If
        `total_units` is zero or there are no events, all labels are `None`.
    """

    if not total_units or not events:
        return [None] * n_items

    labels: List[Optional[Hashable]] = []
    for idx in range(n_items):
        pos = idx * total_units / n_items
        label = None
        for event_pos, event_label in events:
            if event_pos <= pos:
                label = event_label
        labels.append(label)

    return labels


def proportional_positions(
    n_items: int, total_units: float, unit_positions: Iterable[float]
) -> List[int]:
    """
    Maps positions given in another unit to item indices, by proportion.

    Each position `p` becomes the index `round(p * n_items / total_units)`,
    with Python's `round` (which rounds halves to the even number). The
    index is the item before which the position falls: 0 is the start, and
    `n_items` the end. It is used for placing lacunae recorded in another
    unit (see `coverage`).

    Example
    ********

    .. code-block:: python

        >>> seqsim.tradition.proportional_positions(10, 4, [0, 1, 4])
        [0, 2, 10]

    :param n_items: The number of items of the witness.
    :param total_units: The size of the witness in the other unit.
    :param unit_positions: The positions in the other unit.
    :return: A list with the item index of each position; empty if
        `total_units` is zero.
    """

    if not total_units:
        return []

    return [round(pos * n_items / total_units) for pos in unit_positions]


def monotone_labels(
    items: Sequence[Hashable],
    label_order: Sequence[Hashable],
    reference: Mapping[Hashable, Hashable],
    prior: Sequence[Optional[Hashable]],
    prior_weight: float = 0.25,
) -> List[Hashable]:
    """
    Labels a sequence of items monotonically, agreeing with a reference.

    Every item gets one label from `label_order`, labels never go backwards
    in that order, and the labelling maximizes the agreement with
    `reference` (the consensus label of each item, scoring 1 per agreement)
    plus `prior_weight` for each agreement with `prior` (for example the
    proportional labels, from `proportional_labels`), which decides the
    items about which the reference says nothing. It is solved by dynamic
    programming over (item, label), in `O(len(items) * len(label_order))`.

    Ties are broken deterministically: the first label of the highest score
    is chosen for the last item, and, going backwards, the earliest label
    among those with the best score so far.

    Example
    ********

    .. code-block:: python

        >>> reference = {"a": "I", "b": "I", "c": "II", "d": "I", "e": "II"}
        >>> seqsim.tradition.monotone_labels("abcde", ["I", "II"], reference, [None] * 5)
        ['I', 'I', 'I', 'I', 'II']

    :param items: The items of the witness, in order.
    :param label_order: The labels, in the order in which they must appear.
    :param reference: A mapping from items to their reference label.
    :param prior: A label (or `None`) for each position of `items`.
    :param prior_weight: The score of an agreement with `prior`. Defaults to
        0.25.
    :return: A list with the label of each item. If there are no labels,
        all items get the empty string.
    """

    n, k = len(items), len(label_order)
    if not n or not k:
        return [label_order[0] if k else ""] * n

    def gain(i: int, j: int) -> float:
        label = label_order[j]
        score = 1.0 if reference.get(items[i]) == label else 0.0
        return score + (prior_weight if prior[i] == label else 0.0)

    best = [[0.0] * k for _ in range(n)]
    back = [[0] * k for _ in range(n)]
    for j in range(k):
        best[0][j] = gain(0, j)
    for i in range(1, n):
        running, arg = best[i - 1][0], 0
        for j in range(k):
            if best[i - 1][j] > running:
                running, arg = best[i - 1][j], j
            best[i][j] = running + gain(i, j)
            back[i][j] = arg

    j = max(range(k), key=lambda c: best[n - 1][c])
    labels: List[Hashable] = [""] * n
    for i in range(n - 1, -1, -1):
        labels[i] = label_order[j]
        j = back[i][j]

    return labels


def reference_labels(
    voters: Iterable[
        Tuple[Sequence[Hashable], float, Sequence[Tuple[float, Hashable]]]
    ],
    rounds: int = 3,
    prior_weight: float = 0.25,
) -> Tuple[Dict[Hashable, Hashable], Dict[Hashable, float]]:
    """
    Finds the consensus label of each item by an iterated majority vote.

    Each voter is a witness that shares the labelling (for example the
    chapter numbering) of the tradition, given as a tuple of its items, its
    size in another unit, and the events giving where each label starts in
    that unit (see `proportional_labels`); voters without events are
    ignored. The labels of each voter are first placed by proportion, and
    each item gets the label most voted for; then, for `rounds` rounds, the
    labels of each voter are refitted to the consensus with
    `monotone_labels` (keeping the proportional labels as prior, and the
    order of the labels in its events), and the vote is repeated.

    Voters are processed in the order given. Ties in a vote go to the label
    that received its first vote earliest, as in `Counter.most_common`.

    Example
    ********

    .. code-block:: python

        >>> voters = [
        ...     ("abcdef", 6, [(0, "I"), (3, "II")]),
        ...     ("abcdef", 6, [(0, "I"), (2, "II")]),
        ...     ("abdcef", 6, [(0, "I"), (3, "II")]),
        ... ]
        >>> label_of, support = seqsim.tradition.reference_labels(voters, rounds=1)
        >>> label_of
        {'a': 'I', 'b': 'I', 'c': 'II', 'd': 'II', 'e': 'II', 'f': 'II'}
        >>> support["c"], support["d"]
        (0.6666666666666666, 0.6666666666666666)
        >>> label_of, support = seqsim.tradition.reference_labels(voters)
        >>> support["c"], support["d"]
        (1.0, 1.0)

    In the first vote, the first voter places "c" in chapter I, and the third
    places "d" in chapter I; refitted to the consensus, all voters agree.

    :param voters: The voters, as `(items, total_units, events)` tuples.
    :param rounds: The number of rounds of the vote. Defaults to 3.
    :param prior_weight: The weight of the proportional labels when
        refitting (see `monotone_labels`). Defaults to 0.25.
    :return: A tuple with a dictionary from each item that received a vote
        to its label, and a dictionary from each such item to the share of
        the votes its label received.
    """

    voters = [
        (list(items), total_units, list(events))
        for items, total_units, events in voters
        if events
    ]
    labels = [
        proportional_labels(len(items), total_units, events)
        for items, total_units, events in voters
    ]

    reference: Dict[Hashable, Hashable] = {}
    support: Dict[Hashable, float] = {}
    for _ in range(rounds):
        votes: Dict[Hashable, Counter] = defaultdict(Counter)
        for (items, _, _), voter_labels in zip(voters, labels):
            for item, label in zip(items, voter_labels):
                if label:
                    votes[item][label] += 1
        reference = {item: vote.most_common(1)[0][0] for item, vote in votes.items()}
        support = {
            item: vote.most_common(1)[0][1] / sum(vote.values())
            for item, vote in votes.items()
        }
        for idx, (items, total_units, events) in enumerate(voters):
            order = [label for _, label in events]
            prior = proportional_labels(len(items), total_units, events)
            labels[idx] = monotone_labels(items, order, reference, prior, prior_weight)

    return reference, support


def fill_labels(
    frame: Sequence[Hashable], label_of: Mapping[Hashable, Hashable]
) -> Dict[Hashable, Hashable]:
    """
    Gives a label to every item of the frame, from its nearest labelled neighbour.

    Items without a label take the label of the nearest preceding labelled
    item in frame order or, before the first labelled item, of the first
    following one. Items of `label_of` keep their label.

    Example
    ********

    .. code-block:: python

        >>> seqsim.tradition.fill_labels("abcde", {"b": "I", "d": "II"})
        {'b': 'I', 'd': 'II', 'a': 'I', 'c': 'I', 'e': 'II'}

    :param frame: The items, in frame order.
    :param label_of: A mapping from (some) items to their label.
    :return: A dictionary with the labels of `label_of`, followed by those
        of the other items of the frame (or of none, if no item of the
        frame has a label).
    """

    filled = dict(label_of)
    last = None
    pending: List[Hashable] = []
    for item in frame:
        if item in label_of:
            last = label_of[item]
            for pending_item in pending:
                filled[pending_item] = last
            pending = []
        elif last is not None:
            filled[item] = last
        else:
            pending.append(item)

    return filled
