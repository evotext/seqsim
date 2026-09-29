# Comparing the contents of manuscripts

Many medieval and early modern books are collections: legendaries,
homiliaries, florilegia, miscellanies, songbooks, anthologies of poems. Their
copies rarely agree on which texts they contain and in which order. Texts are
added and dropped, gatherings are lost or bound in a different order, the
first and last leaves are damaged, and compilers excerpt, reorganize and
combine their sources.

This tutorial compares the contents of six invented copies of a legendary, a
collection of saints' lives, separating three questions: which texts are
shared, how their order differs, and what was lost or added.

## The manuscripts

Each manuscript is a list of the texts it contains, in order. Texts are
identified by the name of the saint; in a real project they would be
identifiers from a catalogue or a repertory (such as the numbers of the
*Bibliotheca Hagiographica Latina*), assigned after identifying each text.

```python
>>> import seqsim
>>> legendary = {
...     "M1": ["Antonius", "Paulus", "Hilarion", "Malchus", "Martinus", "Ambrosius",
...            "Augustinus", "Hieronymus", "Benedictus", "Gregorius", "Nicolaus", "Silvester"],
...     # M2: a copy of M1 that lost a gathering with three lives
...     "M2": ["Antonius", "Paulus", "Hilarion", "Malchus", "Martinus",
...            "Benedictus", "Gregorius", "Nicolaus", "Silvester"],
...     # M3: a copy of M1 with the first leaves lost, and a life added at the end
...     "M3": ["Hilarion", "Malchus", "Martinus", "Ambrosius", "Augustinus", "Hieronymus",
...            "Benedictus", "Gregorius", "Nicolaus", "Silvester", "Remigius"],
...     # M4: a copy of M1 with the first two quires bound in the opposite order
...     "M4": ["Martinus", "Ambrosius", "Augustinus", "Hieronymus", "Antonius", "Paulus",
...            "Hilarion", "Malchus", "Benedictus", "Gregorius", "Nicolaus", "Silvester"],
...     # M5: a copy of M4 that moved the life of Nicholas to the beginning
...     "M5": ["Nicolaus", "Martinus", "Ambrosius", "Augustinus", "Hieronymus", "Antonius",
...            "Paulus", "Hilarion", "Malchus", "Benedictus", "Gregorius", "Silvester"],
...     # M6: an anthology that excerpted the desert fathers and added other texts
...     "M6": ["Antonius", "Paulus", "Hilarion", "Malchus", "Pachomius", "Macarius"],
... }
>>> names = sorted(legendary)
>>> def show(measure, digits=0, **kwargs):
...     print("    " + "".join(f"{name:>6}" for name in names))
...     for x in names:
...         row = [measure(legendary[x], legendary[y], **kwargs) for y in names]
...         print(f"{x:>4}" + "".join(f"{value:>6.{digits}f}" for value in row))
```

## Which texts are shared?

The Jaccard dissimilarity compares the sets of texts, ignoring their order:

```python
>>> show(seqsim.token.jaccard_dissim, digits=2)
        M1    M2    M3    M4    M5    M6
  M1  0.00  0.25  0.23  0.00  0.00  0.71
  M2  0.25  0.00  0.46  0.25  0.25  0.64
  M3  0.23  0.46  0.00  0.23  0.23  0.87
  M4  0.00  0.25  0.23  0.00  0.00  0.71
  M5  0.00  0.25  0.23  0.00  0.00  0.71
  M6  0.71  0.64  0.87  0.71  0.71  0.00
```

M1, M4 and M5 contain exactly the same texts, so their dissimilarity is zero;
the differences between them are only a matter of order. M6 is far from all
the others, but the Jaccard dissimilarity hides the fact that most of its
texts come from them. Containment is directional and answers "how much of
M6 is found in M1?":

```python
>>> seqsim.token.containment(legendary["M6"], legendary["M1"])
0.6666666666666666
>>> seqsim.token.containment(legendary["M1"], legendary["M6"])
0.3333333333333333
```

Two thirds of M6 are found in M1, but only a third of M1 is found in M6. With
`size=2`, pairs of consecutive texts are compared, which also requires them to
be in the same order: M6 copied a block of four consecutive lives.

```python
>>> seqsim.token.containment(legendary["M6"], legendary["M1"], size=2)
0.6
```

## How different is the order?

The `order` module measures how differently the manuscripts order the texts
they share. The Ulam distance counts the texts that must be moved one by one
(plus those added or removed):

```python
>>> show(seqsim.order.ulam_dist)
        M1    M2    M3    M4    M5    M6
  M1     0     3     3     4     5    10
  M2     3     0     6     4     5     7
  M3     3     6     0     5     6    13
  M4     4     4     5     0     1    10
  M5     5     5     6     1     0    10
  M6    10     7    13    10    10     0
```

For M4, the Ulam distance counts four texts moved, but the change was a single
event: two quires bound in the opposite order. The block interchange measure
counts exchanges of two blocks of texts, and recognizes a single event:

```python
>>> seqsim.order.block_interchange_dissim(legendary["M1"], legendary["M4"])
1.0
>>> seqsim.order.block_interchange_dissim(legendary["M1"], legendary["M5"])
2.0
```

The breakpoint dissimilarity counts pairs of consecutive texts that are no
longer consecutive, and is often a more robust measure of how much the order
was disturbed:

```python
>>> seqsim.order.breakpoint_dissim(legendary["M1"], legendary["M4"])
3.0
>>> seqsim.order.breakpoint_dissim(legendary["M1"], legendary["M5"])
5.0
```

When many rearrangements have happened, the breakpoints saturate, as a
rearrangement can break adjacencies that were already broken. The IEBP
estimate, used for the order of the *Canterbury Tales* by Spencer et al.
(2003), estimates the number of rearrangements (moves of one or more texts)
from the breakpoints of the shared texts:

```python
>>> seqsim.order.iebp_estimate(legendary["M1"], legendary["M4"])
1.0
>>> seqsim.order.iebp_estimate(legendary["M1"], legendary["M5"])
2.0
```

## What was lost, and where?

The Levenshtein distance counts a lost gathering of three texts as three
events, and treats the loss of the first leaves like any other loss:

```python
>>> seqsim.edit.levenshtein_dist(legendary["M1"], legendary["M2"])
3.0
>>> seqsim.edit.levenshtein_dist(legendary["M1"], legendary["M3"])
3.0
```

The stemmatological dissimilarity was designed for these cases: a block of up
to `max_del_len` consecutive texts can be lost (or added) as a single event,
and losses at the beginning and end of a manuscript, which are the most
exposed to damage, cost half as much.

```python
>>> seqsim.edit.stemmatological_dissim(legendary["M1"], legendary["M2"])
1.0
>>> seqsim.edit.stemmatological_dissim(legendary["M1"], legendary["M3"])
1.5
```

The loss of the gathering in M2 counts as a single event. For M3, the text
added at the end costs half a unit, but the loss of the first two texts costs
a full unit: with the default fragile regions of 10% of each manuscript, only
the first text of M1 lies in its fragile beginning. The size of the fragile
regions should reflect the material: for manuscripts that often lost their
first quire, a larger region is appropriate.

```python
>>> seqsim.edit.stemmatological_dissim(legendary["M1"], legendary["M3"], frag_start=20.0)
1.0
```

## Putting it together

No single number describes how two collections differ. A useful strategy is
to compute several matrices, each answering one question, and to compare
them. For example, for each pair of manuscripts, the proportion of shared
texts and the number of rearrangements of those texts:

```python
>>> for x, y in [("M1", "M2"), ("M1", "M4"), ("M4", "M5"), ("M1", "M6")]:
...     shared = 1 - seqsim.token.jaccard_dissim(legendary[x], legendary[y])
...     moves = seqsim.order.iebp_estimate(legendary[x], legendary[y])
...     print(f"{x}-{y}: {shared:.0%} shared, {moves:.0f} rearrangement(s)")
M1-M2: 75% shared, 0 rearrangement(s)
M1-M4: 100% shared, 1 rearrangement(s)
M4-M5: 100% shared, 1 rearrangement(s)
M1-M6: 29% shared, 0 rearrangement(s)
```

M2 and M6 agree with M1 in the order of the texts they share: their
differences are losses, additions, and selection. M4 and M5 contain the same
texts as M1, rearranged.
