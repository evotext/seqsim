# Token measures

The `seqsim.token` module compares the elements, or the short sub-sequences,
that two sequences have in common, regardless of where they occur. These
measures answer questions such as "which texts do two manuscripts share?" or
"how much of this collection is found in that one?".

## Jaccard

`token.jaccard_dissim(x, y, *, normal=False)`

One minus the Jaccard index of the **sets** of elements:
$1 - |X \cap Y| / |X \cup Y|$. Order and repetitions are ignored.

```python
>>> import seqsim
>>> ms_a = ["Vita Antonii", "Vita Pauli", "Vita Hilarionis", "Vita Malchi"]
>>> ms_b = ["Vita Malchi", "Vita Pauli", "Vita Martini"]
>>> seqsim.token.jaccard_dissim(ms_a, ms_b)
0.6
```

Properties
: Symmetric and satisfies the triangle inequality (the Jaccard distance is a
  metric on sets), but different sequences with the same elements have a
  dissimilarity of zero.

When to use
: For the overlap of the contents of manuscripts, ignoring their order.

References
: Jaccard (1912); Tan et al. (2005).

(sorensen-dice)=
## Sørensen-Dice

`token.sorensen_dissim(x, y, *, normal=False)`

One minus the Sørensen-Dice coefficient of the **multisets** of elements,
$1 - 2|X \cap Y| / (|X| + |Y|)$, where repeated elements are counted.

```python
>>> seqsim.token.sorensen_dissim(ms_a, ms_b)
0.4285714285714286
```

Properties
: Symmetric, but different sequences can score zero and the triangle
  inequality does not hold.

References
: Sørensen (1948); Dice (1945).

## Sub-sequence Jaccard

`token.subseq_jaccard_dissim(x, y, *, normal=False)`

For each length $n$ from 1 to the length of the longest sequence, the Jaccard
index of the multisets of contiguous sub-sequences of $n$ elements; the
similarity is the mean of these indices weighted by $n$, so that longer
shared sub-sequences count more, and the dissimilarity is one minus it.

```python
>>> seqsim.token.subseq_jaccard_dissim("abc", "bcde")
0.91
```

Properties
: Symmetric and zero only for identical sequences. No counterexample to the
  triangle inequality is known, but it is not proven.

## q-grams

`token.qgram_dissim(x, y, *, q=2, pad=True, normal=False)`

The q-gram distance of Ukkonen (1992): the sum, over all contiguous
sub-sequences of `q` elements, of the difference of their number of
occurrences in each sequence. It is fast (linear in the length of the
sequences) and a lower bound for the edit distance, which makes it useful for
long texts and as a filter. By default, sequences are padded with `q - 1`
boundary markers on each side, so that the first and last elements are
counted as often as the others; `pad=False` gives Ukkonen's original
definition.

```python
>>> seqsim.token.qgram_dissim("01000", "001111", pad=False)
5.0
>>> text_a = "in principio erat verbum et verbum erat apud deum".split()
>>> text_b = "in principio erat sermo et sermo erat apud deum".split()
>>> seqsim.token.qgram_dissim(text_a, text_b, q=2, normal=True)
0.4
```

Properties
: Symmetric and satisfies the triangle inequality, but different sequences
  can have the same q-grams (e.g., `"abaca"` and `"acaba"`).

References
: Ukkonen (1992).

## Tversky

`token.tversky_simil(x, y, *, alpha=0.5, beta=0.5, normal=False)`

The Tversky index generalizes the Jaccard and Sørensen-Dice coefficients by
weighting the elements unique to each sequence differently:
$|X \cap Y| / (|X \cap Y| + \alpha |X - Y| + \beta |Y - X|)$, on multisets.
With `alpha = beta = 0.5` it is the Sørensen-Dice coefficient; with
`alpha = beta = 1`, the Jaccard index. With different weights it is
directional: with `alpha=1` and `beta=0`, it measures how much of `x` is
found in `y`.

```python
>>> seqsim.token.tversky_simil(["T1", "T2"], ["T1", "T2", "T3", "T4"], alpha=1.0, beta=0.0)
1.0
>>> seqsim.token.tversky_simil(["T1", "T2"], ["T1", "T2", "T3", "T4"])
0.6666666666666666
```

References
: Tversky (1977).

## Containment

`token.containment(x, y, *, size=1)`

The proportion of the contiguous sub-sequences ("shingles") of `size`
elements of `x` that are also found in `y`, counted as multisets (Broder,
1997). It is directional by design, answering questions such as "is
manuscript `x` an excerpt of manuscript `y`?". With `size` larger than one,
it also requires the shared elements to appear in the same local order.

```python
>>> excerpt = ["T3", "T4", "T5"]
>>> source = ["T1", "T2", "T3", "T4", "T5", "T6"]
>>> seqsim.token.containment(excerpt, source)
1.0
>>> seqsim.token.containment(source, excerpt)
0.5
>>> seqsim.token.containment(["T5", "T4", "T3"], source, size=2)
0.0
```

References
: Broder (1997).
