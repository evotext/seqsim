# Case study: the *Dietsche Catoen*

This case study applies `seqsim` to a real tradition: the *Dietsche Catoen*,
the thirteenth-century Middle Dutch translation of the *Distichs of Cato*, a
collection of moral sayings widely used in schools. The data used here cover
thirteen manuscripts and six early printed editions, which differ in the
selection and order of the strophes as well as in their wording. We ask two
questions of the witnesses: how do they compare in the **order of the
strophes**, and how do they compare in their **text**?

## The data

The data come from the diplomatic transcriptions published by Moors,
Voorneveld and van Dalen-Oskam (2025), based on the edition by Van Buuren
(1998). From these transcriptions, we derived a small file with, for each
witness, the list of the strophes of the *Distichs* in the order of the
witness, and the text of three strophes. The file,
{download}`dietsche_catoen.json <../data/dietsche_catoen.json>`, is
distributed under the same license as the transcriptions, CC BY-SA 4.0 (see
[the data folder](https://github.com/evotext/seqsim/tree/main/docs/data) for
details). Download it to follow this tutorial.

```python
>>> import json
>>> import seqsim
>>> with open("dietsche_catoen.json", encoding="utf-8") as handle:
...     data = json.load(handle)
>>> orders = data["orders"]
>>> sorted(orders)
['A', 'B', 'C', 'D', 'G', 'H', 'L', 'M', 'Me', 'P', 'R', 'b', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6']
>>> {witness: len(strophes) for witness, strophes in orders.items()}
{'A': 92, 'B': 51, 'C': 80, 'D': 56, 'G': 36, 'H': 60, 'L': 26, 'M': 48, 'Me': 20, 'P': 69, 'R': 16, 'b': 28, 'd1': 36, 'd2': 36, 'd3': 36, 'd4': 36, 'd5': 36, 'd6': 36}
```

Strophes are identified by book and number (`"I,19"` is the nineteenth
strophe of the first book, `"II,Pa"` a part of the prologue of the second
book). The witnesses differ greatly in size: some are complete, others are
fragments (R and Me) or contain a selection. The first strophes of witness
A and of the first printed edition, d1:

```python
>>> orders["A"][:12]
['I,01', 'I,02', 'I,03', 'I,04', 'I,05', 'I,06', 'I,07', 'I,08', 'I,09', 'I,10', 'I,11', 'I,12']
>>> orders["d1"][:16]
['I,01', 'I,02', 'I,03', 'I,04', 'I,05', 'I,08', 'I,12', 'I,15', 'I,16', 'I,18', 'I,21', 'I,22', 'I,19', 'I,23', 'I,26', 'I,27']
```

The printed edition selects strophes, and places I,19 after I,22.

## The order of the strophes

We compare the witnesses in two ways: the proportion of strophes they do not
share (the Jaccard dissimilarity), and how differently they order the
strophes they share. For the latter, the normalized Ulam distance counts the
strophes that must be moved, added, or removed, relative to the size of the
witnesses.

```python
>>> witnesses = sorted(orders)
>>> def show(measure, **kwargs):
...     print("    " + "".join(f"{w:>5}" for w in witnesses))
...     for x in witnesses:
...         row = [measure(orders[x], orders[y], **kwargs) for y in witnesses]
...         print(f"{x:>4}" + "".join(f"{value:5.2f}" for value in row))
>>> show(seqsim.order.ulam_dist, normal=True)
        A    B    C    D    G    H    L    M   Me    P    R    b   d1   d2   d3   d4   d5   d6
   A 0.00 0.77 0.71 0.78 0.80 0.78 0.90 0.77 0.86 0.78 0.90 0.82 0.80 0.80 0.80 0.80 0.80 0.80
   B 0.77 0.00 0.68 0.65 0.76 0.78 0.80 0.70 0.75 0.81 0.90 0.73 0.76 0.76 0.76 0.76 0.76 0.76
   C 0.71 0.68 0.00 0.77 0.78 0.83 0.89 0.78 0.87 0.85 0.90 0.80 0.78 0.78 0.78 0.78 0.78 0.78
   D 0.78 0.65 0.77 0.00 0.74 0.75 0.86 0.70 0.70 0.78 0.85 0.55 0.74 0.74 0.74 0.74 0.74 0.74
   G 0.80 0.76 0.78 0.74 0.00 0.78 0.78 0.72 0.77 0.77 0.56 0.78 0.00 0.00 0.00 0.00 0.00 0.00
   H 0.78 0.78 0.83 0.75 0.78 0.00 0.86 0.72 0.82 0.84 0.90 0.81 0.78 0.78 0.78 0.78 0.78 0.78
   L 0.90 0.80 0.89 0.86 0.78 0.86 0.00 0.81 0.95 0.82 0.98 0.80 0.78 0.78 0.78 0.78 0.78 0.78
   M 0.77 0.70 0.78 0.70 0.72 0.72 0.81 0.00 0.82 0.81 0.86 0.65 0.72 0.72 0.72 0.72 0.72 0.72
  Me 0.86 0.75 0.87 0.70 0.77 0.82 0.95 0.82 0.00 0.88 0.65 0.72 0.77 0.77 0.77 0.77 0.77 0.77
   P 0.78 0.81 0.85 0.78 0.77 0.84 0.82 0.81 0.88 0.00 0.89 0.87 0.77 0.77 0.77 0.77 0.77 0.77
   R 0.90 0.90 0.90 0.85 0.56 0.90 0.98 0.86 0.65 0.89 0.00 0.87 0.56 0.56 0.56 0.56 0.56 0.56
   b 0.82 0.73 0.80 0.55 0.78 0.81 0.80 0.65 0.72 0.87 0.87 0.00 0.78 0.78 0.78 0.78 0.78 0.78
  d1 0.80 0.76 0.78 0.74 0.00 0.78 0.78 0.72 0.77 0.77 0.56 0.78 0.00 0.00 0.00 0.00 0.00 0.00
  d2 0.80 0.76 0.78 0.74 0.00 0.78 0.78 0.72 0.77 0.77 0.56 0.78 0.00 0.00 0.00 0.00 0.00 0.00
  d3 0.80 0.76 0.78 0.74 0.00 0.78 0.78 0.72 0.77 0.77 0.56 0.78 0.00 0.00 0.00 0.00 0.00 0.00
  d4 0.80 0.76 0.78 0.74 0.00 0.78 0.78 0.72 0.77 0.77 0.56 0.78 0.00 0.00 0.00 0.00 0.00 0.00
  d5 0.80 0.76 0.78 0.74 0.00 0.78 0.78 0.72 0.77 0.77 0.56 0.78 0.00 0.00 0.00 0.00 0.00 0.00
  d6 0.80 0.76 0.78 0.74 0.00 0.78 0.78 0.72 0.77 0.77 0.56 0.78 0.00 0.00 0.00 0.00 0.00 0.00
```

The most striking result is that the six printed editions share exactly the
same selection and order of strophes, and that one manuscript, **G**, shares
it too. The fragment R is closest to this group. Among the other
manuscripts, the closest pairs are D and b, and D and B.

Most of the distances between manuscripts are high because the witnesses
select different strophes. To separate selection from order, we can compare
the proportion of shared strophes with the number of rearrangements of the
strophes they share, estimated with IEBP:

```python
>>> for x, y in [("A", "C"), ("D", "b"), ("B", "D"), ("G", "d1")]:
...     shared = 1 - seqsim.token.jaccard_dissim(orders[x], orders[y])
...     moves = seqsim.order.iebp_estimate(orders[x], orders[y])
...     print(f"{x}-{y}: {shared:.0%} of the strophes shared, {moves:.0f} rearrangement(s)")
A-C: 43% of the strophes shared, 4 rearrangement(s)
D-b: 45% of the strophes shared, 0 rearrangement(s)
B-D: 49% of the strophes shared, 10 rearrangement(s)
G-d1: 100% of the strophes shared, 0 rearrangement(s)
```

D and b share less than half of their strophes, but in the same order: their
differences are a matter of selection. B and D share about the same
proportion, but in a different order. Such contrasts are hard to see from a
single measure.

## The text

The file also holds the text of three strophes (I,01, I,03, and I,19), which
we compare word by word in the fourteen witnesses that contain all three:

```python
>>> strophes = data["strophes"]
>>> complete = [w for w in witnesses if all(w in strophes[s] for s in strophes)]
>>> complete
['A', 'B', 'C', 'D', 'G', 'H', 'M', 'P', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6']
>>> strophes["I,19"]["A"]
['Hijs sotter dan en kuekenoet', 'Die hoept op anders mānes doet', 'Want alle de liede ionc ende out', 'Sijn haers leuens euē ghewout']
>>> strophes["I,19"]["d1"]
['Hy is sotter dan een queken hoot', 'Die hopet op een anders doot', 'Want beyde die lieden ionck ende out', 'Sijn haers leuens onghewout.']
```

The transcriptions are diplomatic: they keep the spelling and the
abbreviation marks of each witness (such as the macron in *mānes*). As we are
interested in the wording rather than the spelling, we normalize the words:
we remove abbreviation marks and punctuation, and treat *u*/*v* and *i*/*j*/*y*
as the same letter.

```python
>>> import re
>>> import unicodedata
>>> def normalize(word):
...     word = unicodedata.normalize("NFD", word.lower())
...     word = "".join(char for char in word if not unicodedata.combining(char))
...     word = re.sub(r"[^a-z]", "", word)
...     return word.replace("v", "u").replace("j", "i").replace("y", "i")
>>> def words(witness):
...     tokens = [normalize(word) for strophe in strophes.values()
...               for line in strophe[witness] for word in line.split()]
...     return [token for token in tokens if token]
>>> texts = {witness: words(witness) for witness in complete}
>>> texts["A"][:8]
['sone', 'maerc', 'wat', 'hic', 'di', 'bediede', 'ende', 'oec']
```

As the witnesses have different lengths, we use the normalized Levenshtein
distance of Yujian and Bo, a true distance suitable for building trees or
networks, and list the two closest witnesses of each one:

```python
>>> for x in complete:
...     closest = sorted(
...         (seqsim.edit.levenshtein_gld_dist(texts[x], texts[y]), y)
...         for y in complete if y != x
...     )[:2]
...     print(x, "->", ", ".join(f"{y} ({value:.2f})" for value, y in closest))
A -> B (0.51), C (0.52)
B -> D (0.41), M (0.41)
C -> P (0.44), M (0.47)
D -> G (0.35), d1 (0.36)
G -> d1 (0.05), d2 (0.08)
H -> D (0.51), B (0.52)
M -> B (0.41), D (0.41)
P -> B (0.42), C (0.44)
d1 -> G (0.05), d2 (0.09)
d2 -> d3 (0.00), d4 (0.00)
d3 -> d2 (0.00), d4 (0.00)
d4 -> d2 (0.00), d3 (0.00)
d5 -> d2 (0.03), d3 (0.03)
d6 -> d5 (0.03), d2 (0.04)
```

The text confirms what the order of the strophes suggested: manuscript G is
very close to the printed editions, and in particular to the first of them,
d1. After normalization, the later editions d2, d3 and d4 are identical in
these strophes, and D is the manuscript closest to the group of G and the
prints.

## Summary

With a few lines of code, two independent kinds of evidence, the order of
the strophes and their wording, point to the same relationships. Such results
are not a stemma: they measure overall similarity, and must be interpreted
with the tools of textual criticism. They are, however, a quick and
reproducible way to explore a tradition, to formulate hypotheses, and to
check them against the whole of the evidence.

## Credits

The data are derived from Moors, S., Voorneveld, N. and van Dalen-Oskam, K.
(2025), "Witnessing Middle Dutch Textual Traditions. Diplomatic Transcriptions
of *Dietsche Catoen*, *Scolastica*, and *Karel ende Elegast*", *Journal of
Open Humanities Data* 11: 43, <https://doi.org/10.5334/johd.328>, dataset
<https://doi.org/10.5281/zenodo.15064631>, licensed under CC BY-SA 4.0. The
transcriptions are based on Van Buuren, F. (ed.) (1998), *De duytschen
Cathoen. Naar de Antwerpse druk van Henrick Eckert van Homberch. Met als
bijlage de andere redacties van de vroegst bekende Middelnederlandse
vertaling der Dicta Catonis*, Hilversum: Verloren.
