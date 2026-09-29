# Getting started

## Installation

`seqsim` requires Python 3.10 or later and has no dependencies. Install it
from PyPI with:

```bash
pip install seqsim
```

To work on the library itself, see [Contributing](contributing.md).

## Comparing two sequences

Each measure is a function that takes two sequences and returns a number.
The functions are organized in modules by family: `edit` (edit distances and
related measures), `order` (the order of shared items), `alignment`,
`token` (shared elements and sub-sequences), `sequence`, and `compression`.

```python
>>> import seqsim
>>> seqsim.edit.levenshtein_dist("kitten", "sitting")
3.0
>>> seqsim.edit.jaro_dissim("kitten", "sitting")
0.25396825396825395
```

A sequence can be a string, a list, a tuple, or any other Python sequence.
Its elements can be anything that can be used as a dictionary key: characters,
words, numbers, tuples, or a mix of them. Two elements are the same if they
are equal (`==`).

```python
>>> seqsim.edit.levenshtein_dist(["the", "black", "cat"], ["the", "cat"])
1.0
>>> seqsim.edit.levenshtein_dist([1, 2, 3], (1, 2, 4))
1.0
```

Note that comparing the strings `"the black cat"` and `"the cat"` compares
characters, while comparing the lists of words compares words. Deciding what
the elements are (characters, words, normalized forms, lemmata, texts) is the
first and most important choice when using the library.

## Normalized values

Many measures return values that grow with the length of the sequences (for
example, the number of edits). Pass `normal=True` to obtain a value in the
range 0 to 1, which allows comparing pairs of sequences of different lengths.

```python
>>> seqsim.edit.levenshtein_dist("kitten", "sitting", normal=True)
0.42857142857142855
```

Measures whose results are always in that range accept `normal` too, so that
all functions can be called in the same way. Parameters other than the two
sequences must always be passed by name:

```python
>>> seqsim.edit.bulk_delete_dist("abcdeXXXXXfghij", "abcdefghij", max_del_len=5)
1.0
```

## The `distance()` wrapper

The `seqsim.distance()` function gives access to all measures of distance and
dissimilarity by name, and accepts more than two sequences, returning the
mean of all pairwise comparisons. Extra keyword arguments are passed to the
measure.

```python
>>> seqsim.distance(["kitten", "sitting"])
3.0
>>> seqsim.distance(["kitten", "sitting", "mitten"], "levenshtein", normal=True)
0.3412698412698412
>>> seqsim.distance(["abcdeXXXXXfghij", "abcdefghij"], "bulk_delete", max_del_len=2)
3.0
```

The available names are the keys of `seqsim.METHODS`:

```python
>>> "ulam" in seqsim.METHODS
True
>>> seqsim.METHODS["ulam"]
<function ulam_dist at ...>
```

## Distance matrices

Most analyses compare every pair of sequences in a collection, producing a
distance matrix that can be passed to clustering, tree-building, or network
methods in other software. A matrix is a short comprehension away:

```python
>>> witnesses = {
...     "A": "in principio erat verbum".split(),
...     "B": "in principio erat sermo".split(),
...     "C": "in principio fuit verbum".split(),
... }
>>> names = sorted(witnesses)
>>> matrix = [
...     [seqsim.edit.levenshtein_dist(witnesses[x], witnesses[y]) for y in names]
...     for x in names
... ]
>>> matrix
[[0.0, 1.0, 1.0], [1.0, 0.0, 2.0], [1.0, 2.0, 0.0]]
```

The [tutorials](tutorials/witnesses.md) show complete examples.

## What next

- [Concepts](concepts.md) explains the difference between distances,
  dissimilarities and similarities, and why it matters.
- [Choosing a method](choosing.md) helps selecting the right measure for a
  research question.
- The [method pages](methods/index.md) describe every measure in detail.
