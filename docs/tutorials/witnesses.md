# Comparing witnesses of a text

This tutorial compares six witnesses of a short text, word by word, to find
which witnesses are closest to each other. It shows how the choice of
elements (raw or normalized words) and of the measure changes the results,
and how to prepare a distance matrix for software that builds trees and
networks.

## The witnesses

The text is the beginning of Genesis in the Latin Vulgate. The six witnesses
are invented, but show the kinds of variation found in real traditions:

- **A** has the standard text;
- **B** has the same text with different spelling (*creauit*, *celum*,
  *uacua*, *tenebre*);
- **C** shares the spelling of B, and omits the words between the two
  occurrences of *super*, a typical eye-skip (*saut du même au même*);
- **D** transposes *deus creavit* and reads *vana* for *vacua*;
- **E** shares the readings of D, omits one *et*, and reads *domini* for *dei*;
- **F** shares the omission of C, and transposes *erat autem*.

```python
>>> import seqsim
>>> witnesses = {
...     "A": "in principio creavit deus caelum et terram terra autem erat inanis "
...          "et vacua et tenebrae super faciem abyssi et spiritus dei ferebatur "
...          "super aquas",
...     "B": "in principio creauit deus celum et terram terra autem erat inanis "
...          "et uacua et tenebre super faciem abyssi et spiritus dei ferebatur "
...          "super aquas",
...     "C": "in principio creauit deus celum et terram terra autem erat inanis "
...          "et uacua et tenebre super aquas",
...     "D": "in principio deus creavit caelum et terram terra autem erat inanis "
...          "et vana et tenebrae super faciem abyssi et spiritus dei ferebatur "
...          "super aquas",
...     "E": "in principio deus creavit caelum et terram terra autem erat inanis "
...          "et vana tenebrae super faciem abyssi et spiritus domini ferebatur "
...          "super aquas",
...     "F": "in principio creauit deus celum et terram terra erat autem inanis "
...          "et uacua et tenebre super aquas",
... }
>>> tokens = {name: text.split() for name, text in witnesses.items()}
>>> names = sorted(tokens)
```

## A first distance matrix

The Levenshtein distance counts the words that must be substituted, added, or
omitted to turn one witness into another. A small helper prints the matrix:

```python
>>> def show(matrix, names, digits=0):
...     print("   " + "".join(f"{name:>7}" for name in names))
...     for name, row in zip(names, matrix):
...         print(f"{name:>3}" + "".join(f"{value:>7.{digits}f}" for value in row))
>>> def distance_matrix(data, measure, **kwargs):
...     return [[measure(data[x], data[y], **kwargs) for y in names] for x in names]
>>> show(distance_matrix(tokens, seqsim.edit.levenshtein_dist), names)
         A      B      C      D      E      F
  A      0      4     11      3      5     13
  B      4      0      7      5      7      9
  C     11      7      0     12     11      2
  D      3      5     12      0      2     14
  E      5      7     11      2      0     13
  F     13      9      2     14     13      0
```

Two groups emerge, A-D-E and B-C-F. However, B differs from A only in
spelling, and yet its distance to A (4) is almost as large as that of E (5),
which has substantive variants: spelling differences count as much as any
other variant.

## Choosing the elements: normalizing spelling

If the research question concerns substantive variants, the spelling should
be normalized before the comparison. `seqsim` compares elements as they are,
so normalization is a preprocessing step:

```python
>>> def normalize(word):
...     return word.replace("v", "u").replace("ae", "e")
>>> normalized = {name: [normalize(word) for word in words] for name, words in tokens.items()}
>>> show(distance_matrix(normalized, seqsim.edit.levenshtein_dist), names)
         A      B      C      D      E      F
  A      0      0      7      3      5      9
  B      0      0      7      3      5      9
  C      7      7      0     10     10      2
  D      3      3     10      0      2     12
  E      5      5     10      2      0     12
  F      9      9      2     12     12      0
```

A and B are now identical, as they have the same substantive text.

## Choosing the measure

The matrix above still has two problems for a textual scholar:

1. The eye-skip in C is a single scribal event, but it counts as seven
   omitted words, making C look far from everyone except F.
2. The transposition *deus creavit* in D counts as two substitutions.

The Damerau-Levenshtein distance counts the transposition of two adjacent
words as a single change:

```python
>>> seqsim.edit.levenshtein_dist(normalized["A"], normalized["D"])
3.0
>>> seqsim.edit.damerau_dist(normalized["A"], normalized["D"])
2.0
```

The global alignment of `seqsim.alignment` can model the eye-skip as a single
event: with a cost for opening a gap, a long omission costs little more than
a short one.

```python
>>> seqsim.alignment.nw_dissim(normalized["A"], normalized["C"])
7.0
>>> seqsim.alignment.nw_dissim(normalized["A"], normalized["C"], gap_open=3.0, gap_extend=0.25)
4.75
>>> seqsim.alignment.nw_dissim(normalized["A"], normalized["D"], gap_open=3.0, gap_extend=0.25)
3.0
```

The omission in C now costs 4.75 instead of 7, while the three substitutions
separating A and D keep their cost.

Which costs are appropriate depends on the tradition and the question; a
useful practice is to check that the conclusions do not change with
reasonable variations of the costs. Note that affine gaps do not satisfy the
triangle inequality (see [Concepts](../concepts.md)).

## Normalized distances for trees and networks

The witnesses have different lengths (C and F are much shorter), so raw
counts are not directly comparable. For a matrix that will be used to build a
tree or network, use a normalized measure that is also a true distance, such
as the normalization of Yujian and Bo:

```python
>>> gld = distance_matrix(normalized, seqsim.edit.levenshtein_gld_dist)
>>> show(gld, names, digits=3)
         A      B      C      D      E      F
  A  0.000  0.000  0.292  0.118  0.192  0.360
  B  0.000  0.000  0.292  0.118  0.192  0.360
  C  0.292  0.292  0.000  0.392  0.400  0.111
  D  0.118  0.118  0.392  0.000  0.082  0.453
  E  0.192  0.192  0.400  0.082  0.000  0.462
  F  0.360  0.360  0.111  0.453  0.462  0.000
```

To find the closest witness to each one:

```python
>>> for x in names:
...     closest = min((y for y in names if y != x), key=lambda y: gld[names.index(x)][names.index(y)])
...     print(x, "->", closest)
A -> B
B -> A
C -> F
D -> E
E -> D
F -> C
```

## Exporting the matrix

Programs for building trees and networks, such as SplitsTree, PHYLIP, or the
`ape` and `phangorn` packages for R, read distance matrices in the PHYLIP
format: the number of taxa, followed by one line per taxon with its name and
distances.

```python
>>> def to_phylip(matrix, names):
...     lines = [str(len(names))]
...     for name, row in zip(names, matrix):
...         lines.append(f"{name:<10}" + " ".join(f"{value:.6f}" for value in row))
...     return "\n".join(lines)
>>> print(to_phylip(gld, names))
6
A         0.000000 0.000000 0.291667 0.117647 0.192308 0.360000
B         0.000000 0.000000 0.291667 0.117647 0.192308 0.360000
C         0.291667 0.291667 0.000000 0.392157 0.400000 0.111111
D         0.117647 0.117647 0.392157 0.000000 0.081633 0.452830
E         0.192308 0.192308 0.400000 0.081633 0.000000 0.461538
F         0.360000 0.360000 0.111111 0.452830 0.461538 0.000000
```

Save the result to a file (for example with
`open("witnesses.phy", "w").write(...)`) and open it in the software of your
choice. In Python, the matrix can also be used with the clustering functions
of SciPy, for example:

```python
from scipy.cluster.hierarchy import average, dendrogram
from scipy.spatial.distance import squareform

tree = average(squareform(gld))
dendrogram(tree, labels=names)
```

Remember that a distance matrix is not a stemma: methods based on distances
group witnesses by overall similarity, while stemmatic reasoning relies on
shared errors. Distances are, however, an excellent tool for exploring a
tradition, checking hypotheses, and dealing with traditions too large for
manual analysis.
