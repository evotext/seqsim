"""
Binary characters for phylogenetic inference: content and adjacency.

Both kinds of character are coded `1`, `0`, or `None` (missing) per witness:

  * **content**: whether a witness has an item, with the coverage model
    deciding where an absence is evidence (see `coverage`);
  * **adjacency**: whether item `x` is immediately followed by item `y` in the
    witness, over the items of a reference frame: `1` if it is, `0` if the
    witness has both but not in that succession, missing if it lacks either.
    This is the binary encoding of gene order used in genome rearrangement
    phylogenetics. An omission that brings `x` and `y` together makes them
    adjacent, so a shared omission is a shared state.

Only characters that vary among the witnesses coding them are kept, and only
if at least two witnesses share each state (`informative`): constant and
singleton characters say nothing about grouping, and the ascertainment
correction of phylogenetic software (such as IQ-TREE) requires their removal.
"""

# Import Python standard libraries
from dataclasses import dataclass
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
from ._coverage import Cell
from ._frame import unique


@dataclass(frozen=True)
class CharacterMatrix:
    """
    A matrix of binary characters: taxa (rows) by characters (columns).

    Rows and columns keep the order in which they were given. Cells are `1`,
    `0`, or `None` (missing). The matrix converts easily to other
    structures, for example to a `pandas` data frame with
    `pandas.DataFrame(matrix.as_dict(), index=matrix.taxa, dtype=object)`.

    Example
    ********

    .. code-block:: python

        >>> matrix = seqsim.tradition.CharacterMatrix(
        ...     taxa=("W1", "W2"), characters=("a", "b"), cells=((1, 0), (None, 1))
        ... )
        >>> matrix.shape
        (2, 2)
        >>> matrix.as_dict()
        {'a': [1, None], 'b': [0, 1]}
        >>> matrix.row("W2")
        (None, 1)

    :param taxa: The names of the taxa (rows).
    :param characters: The names of the characters (columns).
    :param cells: The rows of the matrix, one tuple of cells per taxon.
    """

    taxa: Tuple[str, ...]
    characters: Tuple[str, ...]
    cells: Tuple[Tuple[Cell, ...], ...]

    def __post_init__(self):
        object.__setattr__(self, "taxa", tuple(self.taxa))
        object.__setattr__(self, "characters", tuple(self.characters))
        object.__setattr__(self, "cells", tuple(tuple(row) for row in self.cells))
        if len(self.cells) != len(self.taxa):
            raise ValueError("The matrix must have one row per taxon.")
        if any(len(row) != len(self.characters) for row in self.cells):
            raise ValueError("Every row must have one cell per character.")

    @classmethod
    def from_columns(
        cls,
        taxa: Sequence[str],
        columns: Mapping[str, Sequence[Cell]],
    ) -> "CharacterMatrix":
        """
        Builds a matrix from a mapping of character names to columns of cells.
        """

        names = list(columns)
        cells = tuple(
            tuple(columns[name][row] for name in names) for row in range(len(taxa))
        )
        return cls(tuple(taxa), tuple(names), cells)

    @property
    def shape(self) -> Tuple[int, int]:
        """
        The number of taxa and of characters.
        """

        return len(self.taxa), len(self.characters)

    def column(self, character: str) -> List[Cell]:
        """
        Returns the cells of a character, in the order of the taxa.
        """

        idx = self.characters.index(character)
        return [row[idx] for row in self.cells]

    def row(self, taxon: str) -> Tuple[Cell, ...]:
        """
        Returns the cells of a taxon, in the order of the characters.
        """

        return self.cells[self.taxa.index(taxon)]

    def as_dict(self) -> Dict[str, List[Cell]]:
        """
        Returns a dictionary from character names to columns, in order.
        """

        return {
            name: [row[idx] for row in self.cells]
            for idx, name in enumerate(self.characters)
        }


def informative(column: Iterable[Cell], min_each: int = 2) -> bool:
    """
    Tells whether a character is informative for grouping.

    A character is informative if at least `min_each` witnesses are coded 1
    and at least `min_each` are coded 0.

    Example
    ********

    .. code-block:: python

        >>> seqsim.tradition.informative([1, 1, 0, 0, None])
        True
        >>> seqsim.tradition.informative([1, 0, 0, 0])
        False

    :param column: The cells of the character.
    :param min_each: The minimum number of witnesses with each state.
        Defaults to 2.
    :return: Whether the character is informative.
    """

    values = list(column)
    return values.count(1) >= min_each and values.count(0) >= min_each


def content_characters(
    coverages: Mapping[str, Mapping[Hashable, Cell]], min_each: int = 2
) -> CharacterMatrix:
    """
    Builds the content characters: which witness has which item.

    The rows are the witnesses, in the order of `coverages`; the columns are
    the items, in the order in which they first appear in the coverages
    (the frame order, when all coverages share the frame), keeping only the
    informative ones (see `informative`). Items are named with `str()`. An
    item missing from the coverage of a witness is coded as missing.

    Example
    ********

    .. code-block:: python

        >>> coverages = {
        ...     "W1": {"a": 1, "b": 1, "c": 1},
        ...     "W2": {"a": 1, "b": 0, "c": 1},
        ...     "W3": {"a": 1, "b": 1, "c": 0},
        ...     "W4": {"a": 1, "b": 0, "c": None},
        ... }
        >>> matrix = seqsim.tradition.content_characters(coverages)
        >>> matrix.as_dict()
        {'b': [1, 0, 1, 0]}

    :param coverages: A mapping from witness names to their coverage, as
        returned by `coverage` or `merge_coverage`.
    :param min_each: The minimum number of witnesses with each state.
        Defaults to 2.
    :return: The content characters, as a `CharacterMatrix`.
    """

    taxa = list(coverages)
    items = unique(item for cells in coverages.values() for item in cells)
    columns = {}
    for item in items:
        cells = [coverages[taxon].get(item) for taxon in taxa]
        if informative(cells, min_each):
            columns[str(item)] = cells

    return CharacterMatrix.from_columns(taxa, columns)


def adjacencies(
    sequence: Sequence[Hashable], universe: Iterable[Hashable]
) -> Set[Tuple[Hashable, Hashable]]:
    """
    Returns the immediate successions of a witness, over the items of a universe.

    The witness is reduced to the first occurrence of each item, and then to
    the items of `universe`, before collecting its pairs of consecutive
    items; items outside the universe thus do not break an adjacency.

    Example
    ********

    .. code-block:: python

        >>> sorted(seqsim.tradition.adjacencies("abXcab", "abc"))
        [('a', 'b'), ('b', 'c')]

    :param sequence: The items of the witness.
    :param universe: The items considered (for example, the frame).
    :return: The set of `(x, y)` adjacencies.
    """

    universe = universe if isinstance(universe, (set, frozenset)) else set(universe)
    order = [item for item in unique(sequence) if item in universe]

    return set(zip(order, order[1:]))


def adjacency_characters(
    witnesses: Mapping[str, Sequence[Hashable]],
    universe: Iterable[Hashable],
    min_each: int = 2,
    parts_of: Optional[Mapping[str, Sequence[Sequence[Hashable]]]] = None,
    separator: str = ">",
    key: Optional[Callable[[Hashable], object]] = None,
) -> CharacterMatrix:
    """
    Builds the adjacency characters: which witness has which succession.

    Each pair `(x, y)` of items of the universe adjacent in at least
    `min_each` witnesses (see `adjacencies`) is a candidate character, coded
    1 for the witnesses where `y` immediately follows `x`, 0 for those that
    have both items but not in that succession, and missing for the others;
    only the informative characters are kept (see `informative`). The
    characters are named `f"{x}{separator}{y}"`, with `str()` applied to the
    items.

    The columns are in the sorted order of the `(x, y)` pairs, or of the
    pairs `(key(x), key(y))` if `key` is given (for items that cannot be
    compared with each other); the order does not depend on hashing, so
    results are reproducible across runs. The order matters: some
    phylogenetic methods, such as those resampling characters, give
    different results for different column orders.

    For a witness made of several parts (`parts_of`, for example a
    manuscript merged with `merge_parts`), successions are read within each
    part: 1 if the items are adjacent in some part, 0 only if both items are
    in the same part but not adjacent there, and missing otherwise, since
    concatenating interleaved parts would invent breaks in order.

    Example
    ********

    .. code-block:: python

        >>> witnesses = {"W1": "abcd", "W2": "abcd", "W3": "acbd", "W4": "acbd"}
        >>> matrix = seqsim.tradition.adjacency_characters(witnesses, "abcd")
        >>> matrix.characters
        ('a>b', 'a>c', 'b>c', 'b>d', 'c>b', 'c>d')
        >>> matrix.column("a>b")
        [1, 1, 0, 0]

    :param witnesses: A mapping from witness names to their items.
    :param universe: The items considered, usually the reference frame.
    :param min_each: The minimum number of witnesses with each state, and
        with the adjacency. Defaults to 2.
    :param parts_of: An optional mapping from witness names to the item
        sequences of their parts, used instead of the items of the witness.
    :param separator: The separator of the two items in character names.
        Defaults to `">"`.
    :param key: An optional function of an item, returning a value used for
        sorting the characters.
    :return: The adjacency characters, as a `CharacterMatrix`.
    """

    parts_of = parts_of or {}
    universe_set = set(universe)
    names = list(witnesses)
    parts = {name: list(parts_of.get(name, [witnesses[name]])) for name in names}
    contents = {name: [set(part) for part in parts[name]] for name in names}
    per_witness = {
        name: set().union(*(adjacencies(part, universe_set) for part in parts[name]))
        for name in names
    }

    candidates: Dict[Tuple[Hashable, Hashable], int] = {}
    for adj in per_witness.values():
        for pair in adj:
            candidates[pair] = candidates.get(pair, 0) + 1

    # Sorted: iterating a set of pairs would order the columns by hash, which
    # for strings changes from run to run (PYTHONHASHSEED)
    if key is None:
        pairs = sorted(candidates)
    else:
        pairs = sorted(candidates, key=lambda pair: (key(pair[0]), key(pair[1])))

    columns = {}
    for x, y in pairs:
        if candidates[(x, y)] < min_each:
            continue
        cells: List[Cell] = []
        for name in names:
            if (x, y) in per_witness[name]:
                cells.append(1)
            elif any(x in content and y in content for content in contents[name]):
                cells.append(0)
            else:
                cells.append(None)
        if informative(cells, min_each):
            columns[f"{x}{separator}{y}"] = cells

    return CharacterMatrix.from_columns(names, columns)


def restrict(
    matrix: CharacterMatrix, taxa: Sequence[str], min_each: int = 2
) -> CharacterMatrix:
    """
    Restricts a matrix to some taxa, keeping the characters still informative.

    Example
    ********

    .. code-block:: python

        >>> matrix = seqsim.tradition.CharacterMatrix(
        ...     ("W1", "W2", "W3", "W4", "W5"),
        ...     ("a", "b"),
        ...     ((1, 1), (1, 0), (0, 1), (0, 0), (0, 1)),
        ... )
        >>> seqsim.tradition.restrict(matrix, ["W4", "W3", "W2", "W1"]).characters
        ('a', 'b')
        >>> seqsim.tradition.restrict(matrix, ["W1", "W2", "W3", "W5"]).characters
        ('a',)

    :param matrix: The character matrix.
    :param taxa: The taxa to keep, in the order of the result.
    :param min_each: The minimum number of taxa with each state. Defaults to 2.
    :return: The restricted matrix.
    """

    rows = [matrix.row(taxon) for taxon in taxa]
    keep = [
        idx
        for idx in range(len(matrix.characters))
        if informative((row[idx] for row in rows), min_each)
    ]

    return CharacterMatrix(
        tuple(taxa),
        tuple(matrix.characters[idx] for idx in keep),
        tuple(tuple(row[idx] for idx in keep) for row in rows),
    )


def concat_characters(matrices: Sequence[CharacterMatrix]) -> CharacterMatrix:
    """
    Joins the characters of several matrices with the same taxa, in order.

    Example
    ********

    .. code-block:: python

        >>> m1 = seqsim.tradition.CharacterMatrix(("W1", "W2"), ("a",), ((1,), (0,)))
        >>> m2 = seqsim.tradition.CharacterMatrix(("W1", "W2"), ("a>b",), ((0,), (None,)))
        >>> seqsim.tradition.concat_characters([m1, m2]).as_dict()
        {'a': [1, 0], 'a>b': [0, None]}

    :param matrices: The matrices, all with the same taxa in the same order.
    :return: A matrix with the characters of all matrices.
    """

    if not matrices:
        return CharacterMatrix((), (), ())

    taxa = matrices[0].taxa
    if any(matrix.taxa != taxa for matrix in matrices):
        raise ValueError("All matrices must have the same taxa, in the same order.")

    return CharacterMatrix(
        taxa,
        tuple(name for matrix in matrices for name in matrix.characters),
        tuple(
            tuple(cell for matrix in matrices for cell in matrix.cells[idx])
            for idx in range(len(taxa))
        ),
    )


def character_blocks(
    names: Sequence[str],
    position: Optional[Mapping[Hashable, int]] = None,
    size: Optional[int] = None,
    label_of: Optional[Mapping[Hashable, Hashable]] = None,
    separator: Optional[str] = ">",
) -> List[Hashable]:
    """
    Assigns each character to a block, from the first item it refers to.

    Characters referring to nearby items are not independent, and methods
    resampling characters (such as a block bootstrap) should resample them
    together. The first item of a character is its name up to the first
    `separator` (the whole name for content characters, or with
    `separator=None`). With `size`, the block is the position of that item
    in `position` (usually the frame order) divided by `size`, rounded down;
    with `label_of`, it is the label of the item. The keys of `position` and
    `label_of` are converted with `str()`, as in character names.

    Example
    ********

    .. code-block:: python

        >>> position = {item: pos for pos, item in enumerate("abcdef")}
        >>> names = ["a", "c>d", "e>a", "f"]
        >>> seqsim.tradition.character_blocks(names, position, size=2)
        [0, 1, 2, 2]
        >>> label_of = dict(zip("abcdef", "IIIIII"[:3] + "JJJ"))
        >>> seqsim.tradition.character_blocks(names, label_of=label_of)
        ['I', 'I', 'J', 'J']

    :param names: The names of the characters.
    :param position: A mapping from items to their position, used with
        `size`.
    :param size: The number of consecutive positions in a block.
    :param label_of: A mapping from items to their label, used instead of
        `position` and `size`.
    :param separator: The separator of items in character names, or `None`
        for names that are single items. Defaults to `">"`.
    :return: A list with the block of each character.
    """

    if (size is None) == (label_of is None):
        raise ValueError("Pass either `size` (with `position`) or `label_of`.")

    first = [name.split(separator)[0] if separator else name for name in names]
    if label_of is not None:
        labels = {str(item): label for item, label in label_of.items()}
        return [labels[item] for item in first]

    if position is None:
        raise ValueError("`size` requires `position`.")
    positions = {str(item): pos for item, pos in position.items()}

    return [positions[item] // size for item in first]
