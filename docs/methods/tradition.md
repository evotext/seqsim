# Traditions

The `seqsim.tradition` subpackage analyses a **tradition**: a set of
witnesses, each an ordered sequence of hashable items (texts, sayings,
strophes, chapters), that share a reference frame. Witnesses of such
traditions differ in which items they have and in what order, as in
collections of sayings, anthologies, or legal and liturgical compilations.
The subpackage prepares the data for phylogenetic software, coding content
and order as binary characters, and deciding where the absence of an item is
evidence. It uses only the Python standard library and does not build trees.

A witness is any sequence of hashable items, a collection is a mapping from
witness names to witnesses, and a character cell is `1`, `0`, or `None`
(missing data). Witnesses are reduced to the first occurrence of each item.
All results are deterministic: they do not depend on the order in which
Python iterates sets, which for strings changes from run to run.

The [tradition tutorial](../tutorials/tradition.md) walks through the whole
process on a real tradition. The examples on this page use a small invented
one:

```python
>>> import seqsim
>>> from seqsim import tradition
>>> witnesses = {
...     "W1": ["a", "b", "c", "d", "e", "f", "g", "h"],
...     "W2": ["a", "b", "c", "d", "e", "f", "g", "h"],
...     "W3": ["a", "b", "e", "c", "d", "f", "h"],
...     "W4": ["a", "e", "c", "d", "f", "h"],
...     "W5": ["a", "c", "d", "e", "f", "g"],
...     "W6": ["c", "d", "f", "g"],
... }
```

## Reference frame

`tradition.consensus_order(sequences)`

Merges the orders of several witnesses into one sequence of all their items:
each witness, in the order given, inserts the items not yet placed directly
after the nearest preceding item of that witness already placed (or at the
start). The first witness fixes the order of its items, so the result
depends on the order of the witnesses; list the most representative first.
When the witnesses largely agree, this is a good approximation to a
consensus ranking.

```python
>>> frame = tradition.consensus_order(witnesses.values())
>>> frame
['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
>>> tradition.consensus_order([["a", "c"], ["b", "c", "x"]])
['b', 'a', 'c', 'x']
```

## Labels

Coverage is defined per section of the frame, with a label (a chapter, a
book) for each item. When the labels cannot be read from the items, they can
be inferred from catalogues that give, for each witness, where each section
starts, in another unit than the items (paragraphs, leaves).

`tradition.proportional_labels(n_items, total_units, events)`
: Labels each of `n_items` positions from `(position, label)` events in
  another unit, of which the witness has `total_units`: item `i` is placed
  at `i * total_units / n_items` and gets the label of the last event at or
  before that position.

`tradition.proportional_positions(n_items, total_units, unit_positions)`
: Converts positions in another unit (for example lacunae) to item indices,
  with Python's `round` (halves to even).

`tradition.monotone_labels(items, label_order, reference, prior, prior_weight=0.25)`
: Labels the items of a witness so that labels never go backwards in
  `label_order`, maximizing the agreement with a `reference` label per item
  plus `prior_weight` per agreement with a `prior` labelling (usually the
  proportional one). Solved by dynamic programming.

`tradition.reference_labels(voters, rounds=3, prior_weight=0.25)`
: The consensus label of each item by an iterated majority vote of the
  witnesses sharing the labelling, each given as `(items, total_units,
  events)`: labels are first placed proportionally, and then refitted to the
  consensus with `monotone_labels` for `rounds` rounds. Returns the label of
  each item, and the share of the votes it won.

`tradition.fill_labels(frame, label_of)`
: Labels the items of the frame without a label from their nearest labelled
  neighbour, preceding or else following.

```python
>>> tradition.proportional_labels(8, 4, [(0, "I"), (2, "II")])
['I', 'I', 'I', 'I', 'II', 'II', 'II', 'II']
>>> voters = [
...     (witnesses["W1"], 4, [(0, "I"), (2, "II")]),
...     (witnesses["W3"], 4, [(0, "I"), (2, "II")]),
...     (witnesses["W4"], 3, [(0, "I"), (1.5, "II")]),
... ]
>>> label_of, support = tradition.reference_labels(voters)
>>> label_of
{'a': 'I', 'b': 'I', 'c': 'I', 'd': 'II', 'e': 'I', 'f': 'II', 'g': 'II', 'h': 'II'}
>>> support["d"]
0.6666666666666666
>>> label_of = tradition.fill_labels(frame, label_of)
```

Two of the three voters place "d" in the second chapter, and "e", which they
move before "c", in the first.

## Coverage

`tradition.coverage(sequence, frame, label_of, gaps=(), min_items=2, min_share=0.06)`

An item absent from a witness is evidence about its transmission only if the
witness covers the place where the item belongs. `coverage` codes every item
of the frame as 1 (present), 0 (absent), or `None` (missing). An absent item
is missing if its label is not covered by the witness, if it lies before the
first or after the last item of the witness in frame order (a mutilated
beginning or end), or if it lies in the frame-order gap where a lacuna of the
witness falls. `gaps` are the indices, in the witness reduced to first
occurrences, of the items before which a lacuna falls.

`tradition.covered_labels(...)` gives the labels a witness covers: those of
which it has at least `min_items` items and at least `min_share` of the
items of the frame. The share threshold turns absences into missing data for
excerpt collections, which take a few items from many sections.

```python
>>> tradition.coverage(witnesses["W5"], frame, label_of)
{'a': 1, 'b': 0, 'c': 1, 'd': 1, 'e': 1, 'f': 1, 'g': 1, 'h': None}
>>> tradition.coverage(witnesses["W5"], frame, label_of, gaps=[1])
{'a': 1, 'b': None, 'c': 1, 'd': 1, 'e': 1, 'f': 1, 'g': 1, 'h': None}
>>> tradition.coverage(witnesses["W6"], frame, label_of)
{'a': None, 'b': None, 'c': 1, 'd': 1, 'e': None, 'f': 1, 'g': 1, 'h': None}
```

W5 ends before "h", which is missing; with a lacuna before its second item,
"b" is missing too. W6 has a single item of the first chapter, which it does
not cover.

`tradition.merge_coverage(coverages, frame)`
: Combines the coverage of the parts of a manuscript: present if any part
  has the item, absent if none has it and some part covers its place,
  missing otherwise.

`tradition.density(coverage)`
: The items present divided by the items present or absent (0.0 if none).

`tradition.merge_parts(parts, max_overlap=0.2, compatible=None)`
: Merges the complementary parts of a codex, given as `(label, items)`:
  from the largest part (by number of distinct items, ties keeping the
  input order), each part sharing less than `max_overlap` of its distinct
  items with the merged content so far (and `compatible` with the largest
  part, if given) is merged. Returns the merged parts in label order, and
  the parts left separate (second copies of a section).

```python
>>> tradition.merge_parts([("A", "abcd"), ("B", "efg"), ("C", "abx")])
([('A', 'abcd'), ('B', 'efg')], [('C', 'abx')])
```

## Characters

Two kinds of binary characters describe a tradition:

- **content**: whether a witness has an item, with the coverage deciding
  where an absence is evidence;
- **adjacency**: whether item `x` is immediately followed by item `y`: 1 if
  it is, 0 if the witness has both but not in that succession, missing if it
  lacks either. This is the binary encoding of gene order used in genome
  rearrangement phylogenetics (e.g., Hu et al. 2014). An omission that
  brings `x` and `y` together makes them adjacent, so a shared omission is
  a shared state.

Only characters that vary are kept, and only if at least `min_each`
witnesses (two by default) share each state (`tradition.informative`):
constant and singleton characters say nothing about grouping, and the
ascertainment correction of programs such as IQ-TREE requires their
removal.

`tradition.content_characters(coverages, min_each=2, universe=None)`
: A matrix of witnesses by items, from a mapping of witness names to
  coverages. Columns follow `universe` (usually the frame) if given, and
  otherwise the order in which items first appear in the coverages, which
  is the frame order when all coverages share the frame; items are named
  with `str()`.

`tradition.adjacency_characters(witnesses, universe, min_each=2, parts_of=None, separator=">", key=None)`
: A matrix of witnesses by adjacencies `"x>y"` over the items of `universe`
  (usually the frame), for the pairs adjacent in at least `min_each`
  witnesses. Columns are in the sorted order of the pairs (or of their
  `key`), so that the output is reproducible. For witnesses made of several
  parts (`parts_of`), successions are read within each part.

```python
>>> coverages = {name: tradition.coverage(seq, frame, label_of) for name, seq in witnesses.items()}
>>> content = tradition.content_characters(coverages)
>>> content.characters
('b', 'g')
>>> content.column("g")
[1, 1, 0, 0, 1, 1]
>>> adjacency = tradition.adjacency_characters(witnesses, frame)
>>> adjacency.characters
('d>e', 'd>f', 'e>c', 'e>f', 'f>h')
>>> adjacency.column("e>c")
[0, 0, 1, 1, 0, None]
```

Both return a `CharacterMatrix`, a frozen data class with `taxa`,
`characters`, and `cells` (one tuple per taxon), which keeps the order of
rows and columns. `as_dict()` returns the columns, so that
`pandas.DataFrame(matrix.as_dict(), index=matrix.taxa, dtype=object)` builds
a data frame.

`tradition.restrict(matrix, taxa, min_each=2)`
: Keeps the given taxa (in that order) and the characters still
  informative among them.

`tradition.concat_characters(matrices)`
: Joins the characters of matrices with the same taxa.

`tradition.character_blocks(names, position=None, size=None, label_of=None, separator=">")`
: A block for each character, from its first item: its position divided by
  `size`, or its label. Characters referring to nearby items are not
  independent, and resampling methods such as a block bootstrap should
  resample them together.

```python
>>> position = {item: idx for idx, item in enumerate(frame)}
>>> tradition.character_blocks(adjacency.characters, position, size=4)
[0, 0, 1, 1, 1]
>>> tradition.character_blocks(adjacency.characters, label_of=label_of)
['II', 'II', 'I', 'I', 'II']
```

## Export

`tradition.to_phylip(matrix)`
: Relaxed PHYLIP, with `?` for missing data, for IQ-TREE or RAxML.

`tradition.to_nexus(matrix, charsets=None)`
: A NEXUS data block with the binary (`restriction`) data type of MrBayes,
  and an optional `sets` block with the partitions.

`tradition.charset_ranges(sizes)` and `tradition.partition_nexus(charsets)`
: The ranges of consecutive partitions, and the NEXUS `sets` block read by
  IQ-TREE.

`tradition.taxon_name(witness_id, replacements=...)`
: Makes a witness name safe for these formats: by default, spaces become
  underscores and `+` becomes `-`. Both writers apply it.

`tradition.mrbayes_calibrations(tip_ages)` and `tradition.mrbayes_constraints(groups)`
: Lines of a MrBayes block with the ages of the tips (`fixed` or `uniform`
  calibrations) and constrained groups of taxa with the age of their
  ancestor. They encode no model: priors, clock, and tree models are left to
  the user.

```python
>>> matrix = tradition.concat_characters([content, adjacency])
>>> print(tradition.to_phylip(matrix), end="")
6 7
W1 1110010
W2 1110010
W3 1001101
W4 0001101
W5 011001?
W6 ?1?1???
>>> charsets = tradition.charset_ranges({"content": 2, "adjacency": 5})
>>> print(tradition.partition_nexus(charsets), end="")
#nexus
begin sets;
  charset content = 1-2;
  charset adjacency = 3-7;
end;
>>> tradition.mrbayes_calibrations({"W1": (500, 500), "W6": (300, 400)})
['  calibrate W1 = fixed(500);', '  calibrate W6 = uniform(300,400);']
```

References
: Hu, Fei; Lin, Yu; Tang, Jijun (2014). "MLGO: phylogeny reconstruction and
  ancestral inference from gene-order data". BMC Bioinformatics 15: 354.
  doi:10.1186/s12859-014-0354-6
: Lewis, Paul O. (2001). "A likelihood approach to estimating phylogeny from
  discrete morphological character data". Systematic Biology 50 (6):
  913–925. doi:10.1080/106351501753462876
