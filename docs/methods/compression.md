# Compression

The `seqsim.compression` module measures similarity through compression: if
two sequences are similar, compressing them together takes little more space
than compressing either of them alone. These measures need no alignment and
capture shared material in any position, but they are only meaningful for
sequences long enough to be compressed, in the order of dozens of elements
or more.

Most of them are instances of the **Normalized Compression Distance** (NCD)
of Cilibrasi and Vitányi (2005), where $C(s)$ is the compressed size of $s$:

$$
\mathrm{NCD}(x, y) = \frac{C(xy) - \min(C(x), C(y))}{\max(C(x), C(y))}
$$

## LZMA NCD

`compression.lzma_ncd_dissim(x, y, *, normal=False)`

The NCD using the LZMA compressor of the Python standard library. Each
distinct element is mapped to a fixed-width code, so that sequences of any
elements can be compressed.

```python
>>> import random
>>> import seqsim
>>> rng = random.Random(1)
>>> text = [rng.choice("abcdefgh") for _ in range(300)]
>>> revised = text[:150] + [rng.choice("abcdefgh") for _ in range(150)]
>>> seqsim.compression.lzma_ncd_dissim(text, list(text)) < 0.1
True
>>> seqsim.compression.lzma_ncd_dissim(text, revised) < 0.5
True
```

Properties
: Symmetric, but identical sequences have a small positive dissimilarity, and
  the triangle inequality holds only approximately. The raw value can
  slightly exceed 1.0; `normal=True` clips it to the range 0 to 1.

References
: Cilibrasi and Vitányi (2005).

## Entropy NCD

`compression.entropy_ncd_dissim(x, y, *, normal=False)`

An NCD where the "compressed size" of a sequence is one plus the Shannon
entropy of the distribution of its elements. It only considers how often
each element occurs, not their order.

```python
>>> seqsim.compression.entropy_ncd_dissim("abc", "bcde")
0.21698794996929216
>>> seqsim.compression.entropy_ncd_dissim("ab", "ba")
0.0
```

Properties
: Symmetric and in the range 0 to 1, but different sequences with the same
  frequencies of elements score zero, and the triangle inequality does not
  hold.

## Lempel-Ziv

`compression.lz76_complexity(x)` and `compression.lz76_dissim(x, y, *, normal=False)`

The Lempel-Ziv (1976) complexity of a sequence is the number of components of
its "exhaustive history": the sequence is read from left to right, and each
new component is the shortest piece that cannot be copied from what was
already read. The distance of Otu and Sayood (2003) measures how much the
complexity of each sequence grows when it is appended to the other:

$$
d^*(x, y) = \frac{\max(c(xy) - c(x), c(yx) - c(y))}{\max(c(x), c(y))}
$$

It works directly on the elements, without mapping them to bytes.

```python
>>> seqsim.compression.lz76_complexity("AACGTACCATTG")
7
>>> s, r, q = "AACGTACCATTG", "CTAGGGACTTAT", "ACGGTCACCAA"
>>> seqsim.compression.lz76_dissim(s, q) < seqsim.compression.lz76_dissim(r, q)
True
```

The last example, from Otu and Sayood (2003), shows that `q` is closer to `s`
than to `r`, as `q` shares the patterns `ACG` and `ACC` with `s`.

Properties
: Symmetric, but identical sequences have a small positive dissimilarity,
  different sequences can score zero (e.g., `"aa"` and `"aaa"`), and the
  triangle inequality does not hold.

References
: Lempel and Ziv (1976); Kaspar and Schuster (1987); Otu and Sayood (2003).
