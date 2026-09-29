# Alignment

The `seqsim.alignment` module aligns two sequences, as done in the collation
of witnesses and in bioinformatics, with costs or scores that can be adapted
to the material.

## Global alignment

`alignment.nw_dissim(x, y, *, sub_cost=..., gap_open=0.0, gap_extend=1.0, normal=False)`

The cost of the best global alignment of the two sequences (the
Needleman-Wunsch algorithm, expressed as a cost to minimize), with affine
gaps computed with the algorithm of Gotoh (1982):

- aligning element `a` with element `b` costs `sub_cost(a, b)`, a function
  you can provide (by default, 0 for equal elements and 1 otherwise);
- a gap of `k` consecutive elements (an omission or an addition) costs
  `gap_open + k * gap_extend`.

With the default costs, the result is the Levenshtein distance. Two common
adaptations for textual traditions are making some substitutions cheaper
(orthographic variants, abbreviations, or known scribal confusions), and
setting `gap_open` so that a single long omission, such as an eye-skip
(*saut du même au même*), costs less than the same number of isolated
omissions.

```python
>>> import seqsim
>>> def spelling_cost(a, b):
...     # Words differing only in u/v or i/j are orthographic variants
...     normalize = lambda word: word.replace("v", "u").replace("j", "i")
...     if a == b:
...         return 0.0
...     return 0.2 if normalize(a) == normalize(b) else 1.0
>>> a = "et vidit deus lucem quod esset bona".split()
>>> b = "et uidit deus lucem quod esset bona".split()
>>> seqsim.alignment.nw_dissim(a, b)
1.0
>>> seqsim.alignment.nw_dissim(a, b, sub_cost=spelling_cost)
0.2
```

With `gap_open`, an omission of four consecutive words costs less than three
scattered omissions:

```python
>>> full = "et vidit deus lucem quod esset bona et divisit".split()
>>> skip = "et vidit deus et divisit".split()
>>> scattered = "et deus lucem esset bona divisit".split()
>>> seqsim.alignment.nw_dissim(full, skip, gap_open=2.0)
6.0
>>> seqsim.alignment.nw_dissim(full, scattered, gap_open=2.0)
9.0
```

Properties
: `sub_cost` must be symmetric, non-negative, and zero for equal elements;
  then the result is symmetric and zero only for identical sequences. It is a
  true distance when `sub_cost` is itself a metric and `gap_open` is zero,
  but not in general (affine gaps violate the triangle inequality). With
  `normal=True`, the cost is divided by the cost of aligning each sequence
  entirely against a gap.

References
: Needleman and Wunsch (1970); Gotoh (1982).

## Local alignment

`alignment.sw_simil(x, y, *, score=..., gap_open=0.0, gap_extend=1.0, normal=False)`

The score of the best local alignment (Smith and Waterman, 1981): the pair of
sub-sequences, one from each sequence, that align best. Aligning `a` with
`b` scores `score(a, b)` (by default, 1 for equal elements and -1 otherwise),
and gaps are penalized as in the global alignment. It finds a shared passage
within two otherwise unrelated sequences, such as a quotation or a text
shared by two different compilations.

```python
>>> florilegium = "sicut dicit augustinus in principio erat verbum et verbum erat apud deum".split()
>>> sermon = "fratres carissimi in principio erat verbum et verbum erat apud deum amen".split()
>>> seqsim.alignment.sw_simil(florilegium, sermon)
9.0
>>> seqsim.alignment.sw_simil(florilegium, sermon, normal=True)
0.75
```

With `normal=True`, the score is divided by the highest of the scores of each
sequence aligned with itself; with the default scores, 1.0 indicates
identical sequences.

Properties
: A similarity, symmetric if `score` is symmetric.

References
: Smith and Waterman (1981); Gotoh (1982).

## Monge-Elkan

`alignment.monge_elkan_simil(x, y, *, inner=..., normal=False)`

A similarity for **sequences of sequences**, such as the lists of titles or
incipits of the texts in two manuscripts, where each title is a string. For
each element of one sequence, it takes the highest `inner` similarity to any
element of the other, and averages these values; the result is the mean of
both directions. By default, `inner` is one minus the normalized Levenshtein
distance. The order of the elements is not taken into account.

```python
>>> titles_a = ["Vita sancti Antonii", "Passio sanctae Agnetis"]
>>> titles_b = ["Passio s. Agnetis virginis", "Vita Antonii abbatis"]
>>> round(seqsim.alignment.monge_elkan_simil(titles_a, titles_b), 3)
0.45
```

Any similarity in the range 0 to 1 can be used as `inner`, for example to
compare titles word by word instead of character by character (see the
[tutorial on identifying texts](../tutorials/identifying.md)).

Properties
: A symmetric similarity in the range 0 to 1, equal to 1.0 for identical
  sequences.

References
: Monge and Elkan (1996).
