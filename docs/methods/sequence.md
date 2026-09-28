# Sequence matching

The `seqsim.sequence` module offers the Ratcliff-Obershelp measure, the
algorithm behind Python's `difflib`.

## Ratcliff-Obershelp

`sequence.ratcliff_obershelp_dissim(x, y, *, normal=False)`

The Ratcliff-Obershelp ("gestalt pattern matching") similarity finds the
longest common block of the two sequences, and then recursively the longest
common blocks to its left and to its right. The similarity is twice the
number of matched elements divided by the total number of elements, and the
dissimilarity is one minus it. As the matching depends on the order of the
arguments, `seqsim` computes it in both orders and keeps the highest
similarity.

```python
>>> import seqsim
>>> seqsim.sequence.ratcliff_obershelp_dissim("kitten", "sitting")
0.3846153846153846
>>> seqsim.sequence.ratcliff_obershelp_dissim([1, 2, 3, 4], [2, 4, 3, 5])
0.5
```

Properties
: Symmetric and zero only for identical sequences, but it does not satisfy
  the triangle inequality.

When to use
: A familiar measure, fast in practice, that favours long shared blocks.

References
: Ratcliff and Metzener (1988).
