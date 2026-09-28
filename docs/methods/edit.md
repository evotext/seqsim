# Edit measures

The `seqsim.edit` module collects measures based on the operations
(*edits*) needed to transform one sequence into the other, together with
related measures based on matching elements and blocks. For witnesses of a
text, an edit is a variant reading: a word substituted, added, omitted, or
transposed. For the contents of manuscripts, an edit is a text replaced,
added, lost, or moved.

## Levenshtein

`edit.levenshtein_dist(x, y, *, normal=False)`

The Levenshtein distance is the minimum number of substitutions, insertions,
and deletions of single elements needed to transform `x` into `y`. It is
computed with the Wagner-Fischer algorithm in time proportional to
`len(x) * len(y)`.

```python
>>> import seqsim
>>> seqsim.edit.levenshtein_dist("kitten", "sitting")
3.0
>>> a = "and god said let there be light".split()
>>> b = "and god sayde let ther be lyght".split()
>>> seqsim.edit.levenshtein_dist(a, b)
3.0
```

Properties
: A true distance (metric). With `normal=True`, the distance is divided by
  the length of the longest sequence; the normalized value is not a metric
  (see [Normalized edit distances](#normalized-edit-distances) for
  alternatives).

When to use
: The default choice for counting differences between witnesses of a text,
  at the level of words or characters.

References
: Levenshtein (1966); Wagner and Fischer (1974).

## Damerau-Levenshtein

`edit.damerau_dist(x, y, *, normal=False)`

The (unrestricted) Damerau-Levenshtein distance also counts the transposition
of two adjacent elements as a single edit, computed with the algorithm of
Lowrance and Wagner (1975).

```python
>>> seqsim.edit.levenshtein_dist("in principio erat".split(), "principio in erat".split())
2.0
>>> seqsim.edit.damerau_dist("in principio erat".split(), "principio in erat".split())
1.0
```

Properties
: A true distance (metric).

When to use
: When inversions of adjacent words are frequent and should count as a
  single variant.

References
: Damerau (1964); Lowrance and Wagner (1975).

## Optimal string alignment

`edit.osa_dissim(x, y, *, normal=False)`

The optimal string alignment (OSA), or "restricted Damerau-Levenshtein
distance", counts transpositions like the Damerau-Levenshtein distance, with
the restriction that no part of the sequence is edited more than once. It is
the algorithm often (mis)named "Damerau-Levenshtein" in software, and was
offered under that name in versions of `seqsim` before 0.4.0.

```python
>>> seqsim.edit.osa_dissim("ca", "abc")
3.0
>>> seqsim.edit.damerau_dist("ca", "abc")
2.0
```

Properties
: Symmetric, and zero only for identical sequences, but it does not satisfy
  the triangle inequality: `osa("ca", "abc")` is 3, while going through `"ac"`
  costs 1 + 1.

When to use
: For compatibility with other software; prefer `damerau_dist` otherwise.

## Indel and LCS

`edit.indel_dist(x, y, *, normal=False)` and `edit.lcs_dist(x, y, *, normal=False)`

The **indel distance** allows only insertions and deletions, so that a
substitution counts as two edits (a deletion and an insertion). It is
`len(x) + len(y) - 2 * LCS(x, y)`, where `LCS` is the length of the longest
common subsequence: the longest sequence of elements found in both, in the
same order, but not necessarily contiguous.

The **LCS distance** of Bakkelund (2009) is the proportion of the longest
sequence not covered by the common subsequence,
$1 - \mathrm{LCS}(x, y) / \max(|x|, |y|)$.

```python
>>> seqsim.common.lcs_length("kitten", "sitting")
4
>>> seqsim.edit.indel_dist("kitten", "sitting")
5.0
>>> seqsim.edit.lcs_dist("kitten", "sitting")
0.4285714285714286
```

Properties
: Both are true distances; the LCS distance is in the range 0 to 1.

When to use
: For lists of contents or texts where a replacement is better modelled as
  something lost and something added than as a variant of the same element.

References
: Needleman and Wunsch (1970); Bakkelund (2009).

## Normalized edit distances

`edit.levenshtein_gld_dist`, `edit.damerau_gld_dist`, `edit.indel_gld_dist`,
and `edit.levenshtein_ned_dist`

Dividing an edit distance by the length of the longest sequence gives a value
between 0 and 1, but breaks the triangle inequality (see
[Concepts](../concepts.md#normalization)). Two normalizations that preserve it
are available.

The normalization of **Yujian and Bo (2007)** divides the distance $d$ by the
lengths of both sequences and $d$ itself:

$$
d_{\mathrm{GLD}}(x, y) = \frac{2 \, d(x, y)}{|x| + |y| + d(x, y)}
$$

It is available for the Levenshtein, Damerau-Levenshtein, and indel
distances, for which it is a true distance (it is the "Steinhaus transform"
of a metric whose distance of each sequence to the empty sequence is its
length). It is 1.0 only when one of the sequences is empty.

The **normalized edit distance of Marzal and Vidal (1993)** is the minimum,
over all ways of editing `x` into `y`, of the number of edits divided by the
total number of operations, including the elements kept unchanged. A longer
sequence of operations with more unchanged elements can therefore have a
lower ratio. With unit costs it is a true distance (Fisman et al., 2022), but
it is computed in time proportional to `len(x) * len(y) * (len(x) + len(y))`,
and is slow for long sequences.

```python
>>> seqsim.edit.levenshtein_gld_dist("kitten", "sitting")
0.375
>>> seqsim.edit.levenshtein_ned_dist("kitten", "sitting")
0.42857142857142855
>>> seqsim.edit.levenshtein_ned_dist("ab", "ba")
0.6666666666666666
```

When to use
: For distance matrices of witnesses of different lengths that will be used
  to infer trees or networks.

References
: Yujian and Bo (2007); Marzal and Vidal (1993); Fisman et al. (2022).

## Bulk delete

`edit.bulk_delete_dist(x, y, *, max_del_len=5, normal=False)`

An edit distance where a block of up to `max_del_len` consecutive elements
can be deleted or inserted in a single operation, while substitutions cost
one. It models the loss (or addition) of a group of texts, such as the loss
of a gathering in a manuscript, as a single event. With `max_del_len=1` it is
the Levenshtein distance.

```python
>>> full = ["T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8"]
>>> lost_gathering = ["T1", "T2", "T6", "T7", "T8"]
>>> seqsim.edit.levenshtein_dist(full, lost_gathering)
3.0
>>> seqsim.edit.bulk_delete_dist(full, lost_gathering, max_del_len=3)
1.0
```

Properties
: A true distance (all operations are symmetric and have positive costs).

References
: Göransson et al. (in preparation).

## Fragile ends

`edit.fragile_ends_dissim(x, y, *, frag_start=10.0, frag_end=10.0, normal=False)`

The Levenshtein distance, but deletions and insertions in the first
`frag_start` and last `frag_end` percent of each sequence cost half as much.
It models manuscripts whose first and last leaves, and thus the texts they
contained, are more likely to be lost.

```python
>>> lost_first = ["T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10"]
>>> full = ["T1"] + lost_first
>>> seqsim.edit.fragile_ends_dissim(full, lost_first)
0.5
>>> seqsim.edit.fragile_ends_dissim(full, full[:4] + full[5:])
1.0
```

Properties
: Symmetric and zero only for identical sequences, but it does not satisfy
  the triangle inequality, as the discount depends on the positions in each
  sequence.

References
: Göransson et al. (in preparation).

## Stemmatological

`edit.stemmatological_dissim(x, y, *, frag_start=10.0, frag_end=10.0, max_del_len=5, normal=False)`

Combines the two previous measures: blocks of up to `max_del_len`
consecutive elements can be deleted or inserted as a single operation, and
blocks lying entirely within the fragile regions at the beginning and end of
each sequence cost half as much. It was developed for comparing the contents
of mixed-content miscellanies, where texts are frequently lost in groups and
at the extremities of the codex.

```python
>>> x = ["T%d" % i for i in range(1, 21)]
>>> y = x[2:10] + x[15:]  # lost the first two texts and a gathering of five
>>> seqsim.edit.levenshtein_dist(x, y)
7.0
>>> seqsim.edit.stemmatological_dissim(x, y)
1.5
```

Properties
: Symmetric and zero only for identical sequences, but it does not satisfy
  the triangle inequality.

References
: Göransson et al. (in preparation).

## Block moves

`edit.block_move_dissim(x, y, *, directional=False, normal=False)`

Following Tichy (1984), `y` is built from `x` by copying blocks of `x`, in any
order and possibly more than once, and adding the elements not found in `x`.
The measure is the minimum number of "cuts" needed: the number of pieces
(blocks, plus one addition for each run of new elements) minus one. Both
sequences are framed by start and end markers, so that identical sequences
need a single piece. By default, the highest value of both directions is
returned.

```python
>>> seqsim.edit.block_move_dissim("abcdefgh", "efghabcd")
3.0
>>> seqsim.edit.block_move_dissim("abc", "abcabc", directional=True)
1.0
```

Properties
: Symmetric and zero only for identical sequences, but it does not satisfy
  the triangle inequality.

When to use
: For sequences where blocks are moved or repeated, such as manuscripts with
  reordered gatherings.

References
: Tichy (1984).

## Greedy String Tiling

`edit.gst_dissim(x, y, *, min_match=2, normal=False)`

Greedy String Tiling covers both sequences with "tiles", common blocks of at
least `min_match` elements, taking the longest available ones first. As the
tiles can appear in any order, it measures how much material is shared in
blocks, regardless of their position, and is widely used for detecting text
reuse. The similarity is the proportion of elements covered by tiles,
$2 \cdot \mathrm{coverage} / (|x| + |y|)$, and the dissimilarity is one minus
it.

```python
>>> seqsim.edit.gst_dissim("abcdefgh", "efghabcd")
0.0
>>> seqsim.edit.gst_dissim("abcdefgh", "efghXbcd")
0.125
```

Properties
: Symmetric, but different sequences can score zero (when they share all
  blocks in a different order), and the triangle inequality does not hold.
  Identical sequences always score zero.

When to use
: For detecting shared passages or groups of texts that were moved as
  blocks. Increase `min_match` to ignore short, accidental coincidences.

References
: Wise (1993); Prechelt et al. (2002).

## Jaro and Jaro-Winkler

`edit.jaro_dissim(x, y, *, normal=False)` and `edit.jaro_winkler_dissim(x, y, *, normal=False)`

The Jaro similarity counts the elements matched within a window around their
position and the transpositions among them; the Jaro-Winkler similarity
raises it for sequences sharing a prefix of up to four elements. The
dissimilarities are one minus the similarities. Both were designed for short
strings such as names.

```python
>>> seqsim.edit.jaro_dissim("MARTHA", "MARHTA")
0.05555555555555547
>>> seqsim.edit.jaro_winkler_dissim("MARTHA", "MARHTA")
0.03888888888888886
```

Properties
: Symmetric and zero only for identical sequences, but they do not satisfy
  the triangle inequality.

When to use
: For short sequences such as personal names, place names, or short titles.

References
: Jaro (1989); Winkler (1990).

## MMCWPA

`edit.mmcwpa_dissim(x, y, *, normal=False)`

The Modified Moving Contracting Window Pattern Algorithm repeatedly removes
the longest block shared by both sequences, accumulating the sum of the
squares of the sizes of the removed blocks (the SSNC). The dissimilarity is
$1 - \sqrt{\mathrm{SSNC}} / (|x| + |y|)$, so that long shared blocks weigh
more than the same number of scattered shared elements.

```python
>>> seqsim.edit.mmcwpa_dissim("kitten", "sitting")
0.5134957445894801
```

Properties
: Symmetric and zero only for identical sequences, but it does not satisfy
  the triangle inequality.

References
: Yang et al. (2001); Tresoldi (2016).

## Birnbaum

`edit.birnbaum_simil(x, y, *, normal=False)` and `edit.birnbaum_dissim(x, y, *, normal=False)`

The Birnbaum similarity, proposed for comparing the contents of mixed-content
miscellanies, aligns the sequences into matching blocks and scores each block
of $n$ elements with $n(n+1)/2$, so that longer runs of texts in the same
order weigh more. The normalized similarity divides the score by that of the
longest sequence compared with itself, and the dissimilarity is one minus the
normalized similarity.

```python
>>> seqsim.edit.birnbaum_simil("kitten", "sitting")
7.0
>>> seqsim.edit.birnbaum_dissim("kitten", "sitting")
0.75
```

Properties
: `birnbaum_dissim` is symmetric and zero only for identical sequences, but
  it does not satisfy the triangle inequality.

References
: Birnbaum (2003).
