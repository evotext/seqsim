# Concepts

## Sequences and elements

`seqsim` compares two **sequences**, ordered collections of **elements**.
What counts as an element is a decision of the researcher, and it changes
what is being measured:

| Elements | Example sequence | What a difference means |
|---|---|---|
| characters | `"verbum"` | a spelling difference |
| words (tokens) | `["in", "principio", "erat", "verbum"]` | a variant reading |
| normalized words or lemmata | `["in", "principium", "esse", "verbum"]` | a substantive variant, ignoring spelling |
| texts in a manuscript | `["Vita Antonii", "Vita Pauli", ...]` | a text added, lost, or moved |
| quires, sections, verses... | `[("Q1", 8), ("Q2", 8), ...]` | a structural change |

Elements can be any hashable Python object: strings, numbers, tuples, or a
mix of them. Two elements are the same if they are equal (`==`). The library
never looks inside an element, so normalizing spelling, lemmatizing, or
assigning identifiers to texts should be done before the comparison.

## Distances, dissimilarities, and similarities

A **distance** (or metric) $d$ is a measure with four properties, for all
sequences $x$, $y$ and $z$:

1. non-negativity: $d(x, y) \geq 0$;
2. identity of indiscernibles: $d(x, y) = 0$ if and only if $x = y$;
3. symmetry: $d(x, y) = d(y, x)$;
4. triangle inequality: $d(x, z) \leq d(x, y) + d(y, z)$.

Many useful measures lack one or more of these properties. For example, the
Jaccard dissimilarity compares the sets of elements, ignoring their order, so
`"ab"` and `"ba"` have a dissimilarity of zero although they are different;
and the Jaro dissimilarity does not satisfy the triangle inequality. `seqsim`
states the properties of each measure in the name of its function:

`_dist`
: A true distance: all four properties hold. For measures that count edits,
  the properties hold for the raw values (see [Normalization](#normalization)).

`_dissim`
: A dissimilarity: identical sequences score `0.0` and higher values indicate
  more different sequences, but the four properties are not all guaranteed.
  The documentation of each function states which properties fail, usually
  with an example.

`_simil`
: A similarity: higher values indicate more similar sequences.

A few functions do not follow this scheme because they are directional by
design, such as `token.containment()` (how much of one sequence is found in
another), or because they are estimators rather than measures, such as
`order.iebp_estimate()`.

The name is a promise that is checked: the test suite verifies the properties
of every measure on thousands of generated sequences, and the claims were
checked against the original publications.

### Why the properties matter

The properties are not a mathematical nicety. Methods that build trees or
networks from a distance matrix, such as neighbour joining, minimum evolution,
or split decomposition, assume that the input behaves like a distance. When
the triangle inequality fails, a matrix can describe "shortcuts" that no tree
or network can represent, leading to negative branch lengths or distorted
topologies. Clustering methods are more tolerant, and dissimilarities are
perfectly adequate for exploring a collection, ranking the witnesses closest
to a given one, or detecting outliers.

As a rule of thumb: prefer a `_dist` measure when the matrix will be used to
infer relationships (a stemma, a tree, a network), and choose freely among
all measures for exploration.

## Normalization

Measures that count operations, such as the Levenshtein distance, return
values that grow with the length of the sequences: two long witnesses will
have more differences than two short ones, even if they are equally close.
Passing `normal=True` returns a value between 0 and 1, usually by dividing by
the largest possible value (for the Levenshtein distance, the length of the
longest sequence).

This common normalization does not preserve the triangle inequality: for the
Levenshtein distance, `"ba"` and `"ab"` have a normalized distance of 1.0, but
both have a normalized distance of only 1/3 to `"bab"`.

```python
>>> seqsim.edit.levenshtein_dist("ba", "ab", normal=True)
1.0
>>> seqsim.edit.levenshtein_dist("ba", "bab", normal=True)
0.3333333333333333
```

When a normalized value must also be a true distance, use the measures
designed for that purpose: `edit.levenshtein_gld_dist()`,
`edit.levenshtein_ned_dist()`, `edit.lcs_dist()`, and their relatives (see
[Edit measures](methods/edit.md)).

```python
>>> seqsim.edit.levenshtein_gld_dist("ba", "ab")
0.6666666666666666
>>> seqsim.edit.levenshtein_gld_dist("ba", "bab")
0.3333333333333333
```

Here the distance between `"ba"` and `"ab"` (2/3) is no longer larger than
the sum of their distances to `"bab"` (1/3 each).

## Symmetry

Some classic algorithms give different results depending on which sequence
comes first, because they match elements greedily from left to right (Jaro,
Ratcliff-Obershelp, MMCWPA, Birnbaum, Greedy String Tiling, Tichy's block
moves). `seqsim` computes these in both directions and keeps the result
indicating the greatest similarity, so that all measures available through
`distance()` are symmetric. Where the directional value is meaningful, it is
available with `directional=True`.

## Empty sequences

Empty sequences are allowed everywhere. Two empty sequences are identical and
have a dissimilarity of zero; an empty sequence and a non-empty one have the
largest possible dissimilarity (1.0 for measures in the range 0 to 1). For
edit distances, the value is the cost of inserting the non-empty sequence.

```python
>>> seqsim.edit.jaro_dissim([], [])
0.0
>>> seqsim.edit.jaro_dissim([], ["a"])
1.0
```

## Repeated elements

Measures of the order of items (the `order` module) are classically defined
for permutations, where each item appears exactly once, such as the list of
texts in a manuscript. When an element is repeated (for example, a text
copied twice), `seqsim` distinguishes occurrences: the first occurrence in one
sequence corresponds to the first occurrence in the other, and so on. For
some measures this is only a heuristic, and the documentation of each
function says so.

## More than two sequences

`seqsim.distance()` accepts any number of sequences and returns the mean of all
pairwise comparisons, a simple summary of the diversity of a group of
witnesses. For most analyses, however, compute the full matrix of pairwise
comparisons, as shown in [Getting started](getting-started.md#distance-matrices).

## Computational cost

Most measures take time proportional to the product of the lengths of the two
sequences, which is fast for texts of thousands of words. A few are slower and
are documented as such, notably `edit.levenshtein_ned_dist()` (cubic) and
`edit.gst_dissim()` on long sequences with many repetitions. The
compression-based measures are only meaningful for sequences long enough to
be compressed.
