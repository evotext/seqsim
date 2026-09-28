# Tutorial: from a tradition to phylogenetic characters

This tutorial prepares the order of the strophes of the *Dietsche Catoen*,
the tradition of the [case study](case-study.md), for phylogenetic software
such as IQ-TREE or MrBayes, with the `seqsim.tradition` subpackage. The
witnesses differ in **which** strophes they have and in **what order**, and
both can be coded as binary characters: whether a witness has a strophe
(content), and whether a strophe is immediately followed by another
(adjacency). Along the way, we decide where the absence of a strophe from a
witness is evidence, and where it is only missing data.

`seqsim.tradition` prepares the data; it does not build trees. The files it
writes can be read by any phylogenetic program.

## The data

We use the file {download}`dietsche_catoen.json <../data/dietsche_catoen.json>`
of the case study, with the strophes of each witness in its order.

```python
>>> import json
>>> import seqsim
>>> from seqsim import tradition
>>> with open("dietsche_catoen.json", encoding="utf-8") as handle:
...     orders = json.load(handle)["orders"]
>>> orders["R"]
['I,01', 'I,02', 'I,03', 'I,04', 'I,05', 'I,08', 'III,16', 'II,21', 'II,17', 'II,24', 'II,26', 'II,28', 'II,29', 'III,01', 'III,21', 'III,21bis']
```

Some witnesses repeat strophes; the functions of `seqsim.tradition` use the
first occurrence of each.

```python
>>> len(orders["A"]), len(tradition.unique(orders["A"]))
(92, 80)
```

## A reference frame

A reference frame is a single order of all the strophes of the tradition.
`consensus_order` merges the orders of the witnesses: the first witness
fixes the order of its strophes, and each following one inserts the
strophes not yet placed after their nearest predecessor in that witness. The
result depends on the order of the witnesses, so we start with the largest
ones, whose order is the most complete.

```python
>>> by_size = sorted(orders, key=lambda witness: -len(set(orders[witness])))
>>> by_size[:6]
['A', 'C', 'D', 'H', 'B', 'P']
>>> frame = tradition.consensus_order(orders[witness] for witness in by_size)
>>> len(frame)
137
>>> frame[:8]
['I,01', 'I,02', 'I,03', 'I,04', 'I,05', 'I,06', 'I,07', 'I,08']
```

Each strophe also needs a label, a section of the work where it belongs. In
traditions with chapters recorded in catalogues, `reference_labels` finds
the chapter of each item by a vote of the witnesses (see the
[method page](../methods/tradition.md)). Here the strophe identifiers give
the book directly:

```python
>>> book = {strophe: strophe.split(",")[0] for strophe in frame}
>>> book["III,21bis"]
'III'
```

## Coverage: where absence is evidence

A strophe missing from a witness tells something about its transmission
only if the witness could have had it. `coverage` codes every strophe of the
frame as present (1), absent (0), or missing (`None`): a strophe is missing
if it belongs to a book the witness does not cover, if it lies before the
first or after the last strophe of the witness in frame order, or inside a
recorded lacuna.

```python
>>> cells = tradition.coverage(orders["b"], frame, book)
>>> sorted(tradition.covered_labels(orders["b"], frame, book))
['I', 'II']
>>> cells["I,01"], cells["I,03"], cells["IV,48"]
(None, 1, None)
>>> from collections import Counter
>>> Counter(cells.values())
Counter({None: 72, 0: 37, 1: 28})
```

Witness b begins with I,03: the first two strophes may have been lost with
the beginning of the manuscript, and are missing, not absent. It has no
strophe of books III and IV, which are missing too.

A witness covers a book if it has at least `min_items` of its strophes
(2 by default) and at least `min_share` of them (6% by default). Raising
`min_share` treats witnesses with only a few strophes of each book as
excerpts, for which absences say little. Witness R has 16 strophes spread
over three books:

```python
>>> Counter(tradition.coverage(orders["R"], frame, book).values())
Counter({0: 76, None: 45, 1: 16})
>>> Counter(tradition.coverage(orders["R"], frame, book, min_share=0.2).values())
Counter({None: 121, 1: 16})
```

Lacunae recorded in the witness are passed as `gaps`, the positions in the
witness (reduced to first occurrences) before which a strophe or more were
lost; `proportional_positions` converts positions given in another unit,
such as leaves, into such indices. If, for example, a leaf was lost in R
between II,29 and III,01 (before its fourteenth strophe, at index 13), the
strophes between these two in frame order are missing:

```python
>>> lost = tradition.coverage(orders["R"], frame, book, gaps=[13])
>>> Counter(lost.values())
Counter({None: 66, 0: 55, 1: 16})
>>> tradition.coverage(orders["R"], frame, book)["II,18"], lost["II,18"]
(0, None)
```

The share of the covered strophes that a witness has, its `density`, tells
complete witnesses from selections:

```python
>>> coverages = {witness: tradition.coverage(orders[witness], frame, book) for witness in orders}
>>> round(tradition.density(coverages["A"]), 3), round(tradition.density(coverages["d1"]), 3)
(0.602, 0.387)
```

## Content characters

`content_characters` turns the coverages into a matrix of witnesses by
strophes, keeping only the strophes informative for grouping: those present
in at least two witnesses and absent from at least two others.

```python
>>> content = tradition.content_characters(coverages)
>>> content.shape
(18, 81)
>>> content.characters[:5]
('I,06', 'I,07', 'I,09', 'I,10', 'I,11')
>>> dict(zip(content.taxa, content.column("I,06")))
{'A': 1, 'B': 1, 'C': 1, 'D': 1, 'G': 0, 'H': 0, 'L': None, 'M': 0, 'Me': 1, 'P': 0, 'R': 0, 'b': 1, 'd1': 0, 'd2': 0, 'd3': 0, 'd4': 0, 'd5': 0, 'd6': 0}
```

The result is a `CharacterMatrix`, which keeps the order of the witnesses
and of the characters, and converts easily to other structures: for example,
`pandas.DataFrame(content.as_dict(), index=content.taxa)`.

## Adjacency characters

The order of the strophes is coded with adjacency characters: for each pair
of strophes `x` and `y` adjacent in at least two witnesses, whether `y`
immediately follows `x` (1), whether the witness has both but not in that
succession (0), or whether it lacks either (`None`). Strophes outside the
frame are ignored, and an omission that brings two strophes together makes
them adjacent, so that a shared omission is a shared state.

```python
>>> adjacency = tradition.adjacency_characters(orders, frame)
>>> adjacency.shape
(18, 39)
>>> adjacency.characters[:4]
('I,04>I,05', 'I,04>I,07', 'I,05>I,08', 'I,06>I,07')
>>> dict(zip(adjacency.taxa, adjacency.column("I,05>I,08")))
{'A': 0, 'B': 0, 'C': 0, 'D': 0, 'G': 1, 'H': None, 'L': None, 'M': 1, 'Me': 0, 'P': None, 'R': 1, 'b': 0, 'd1': 1, 'd2': 1, 'd3': 1, 'd4': 1, 'd5': 1, 'd6': 1}
```

The columns are in the sorted order of the pairs of strophes, so that the
output does not change from run to run. For a manuscript made of several
parts (collections bound together, or interleaved sections), pass the parts
with `parts_of`: successions are then read within each part.

## Selecting witnesses

To analyse only the manuscripts, without the printed editions, `restrict`
keeps some witnesses and drops the characters that are no longer
informative among them:

```python
>>> manuscripts = [witness for witness in orders if not witness.startswith("d")]
>>> content_ms = tradition.restrict(content, manuscripts)
>>> adjacency_ms = tradition.restrict(adjacency, manuscripts)
>>> content_ms.shape, adjacency_ms.shape
((12, 80), (12, 30))
```

## Export

The content and adjacency characters are usually analysed together, as two
partitions with their own model. `concat_characters` joins them, and
`charset_ranges` gives the range of each partition:

```python
>>> combined = tradition.concat_characters([content_ms, adjacency_ms])
>>> charsets = tradition.charset_ranges(
...     {"content": content_ms.shape[1], "adjacency": adjacency_ms.shape[1]}
... )
>>> charsets
{'content': (1, 80), 'adjacency': (81, 110)}
```

`to_phylip` writes relaxed PHYLIP, with `?` for missing data, for IQ-TREE or
RAxML, and `partition_nexus` the partition file of IQ-TREE:

```python
>>> phylip = tradition.to_phylip(combined)
>>> print("\n".join(line[:60] for line in phylip.splitlines()[:3]))
12 110
A 1111111111111111000000000111111110111111100011000000110001
B 1111010111111111110111100011011101111000000001110000101101
>>> print(tradition.partition_nexus(charsets), end="")
#nexus
begin sets;
  charset content = 1-80;
  charset adjacency = 81-110;
end;
```

`to_nexus` writes a NEXUS data block for MrBayes, with the binary
("restriction") data type and, optionally, the partitions:

```python
>>> nexus = tradition.to_nexus(combined, charsets)
>>> print("\n".join(line[:60] for line in nexus.splitlines()[:7]))
#NEXUS
begin data;
  dimensions ntax=12 nchar=110;
  format datatype=restriction missing=? gap=-;
  matrix
    A   1111111111111111000000000111111110111111100011000000
    B   1111010111111111110111100011011101111000000001110000
>>> print("\n".join(nexus.splitlines()[-5:]))
end;
begin sets;
  charset content = 1-80;
  charset adjacency = 81-110;
end;
```

Save the text with, for example,
`pathlib.Path("catoen.phy").write_text(phylip, encoding="utf-8")` (and the
partitions as `catoen.nex`), and run
IQ-TREE with a binary model and ascertainment correction (all characters
are variable), such as `iqtree2 -s catoen.phy -p catoen.nex -m GTR2+FO+ASC`.
Taxon names are made safe for these programs with `taxon_name`, which
replaces spaces and `+` signs.

## Comparing pairs of witnesses

The characters describe the whole tradition; the order measures of
`seqsim.order` compare two witnesses directly. `kendall_tau_simil` is
Kendall's correlation between the orders of the strophes two witnesses
share, and `breakpoint_simil`, without boundaries and on the shared
strophes, the share of their successions that is preserved, which is more
sensitive to the move of a few blocks:

```python
>>> seqsim.order.kendall_tau_simil(orders["A"], orders["d1"])
0.86
>>> shared_a, shared_d1 = seqsim.order.restrict_to_shared(orders["A"], orders["d1"])
>>> seqsim.order.breakpoint_simil(shared_a, shared_d1, boundaries=False)
0.7083333333333334
```
