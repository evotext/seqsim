# Order measures

The `seqsim.order` module measures how differently two sequences order the
items they share. The measures come from two fields where the same question
arises: the comparison of rankings (Kendall tau, Spearman footrule), and the
study of genome rearrangements, where genes are reordered on chromosomes by
events that move or exchange blocks of genes (breakpoints, block
interchanges). In textual scholarship, they apply naturally to the order of
texts in manuscripts, of poems in anthologies, of tales in collections, or of
chapters in versions of a work.

These measures are classically defined for **permutations**: two sequences of
the same distinct items. `seqsim` generalizes them to any pair of sequences:

- items present in only one sequence (a text copied in only one manuscript)
  are handled as described for each measure, usually counting as one
  insertion or deletion;
- repeated items (a text copied twice) are distinguished by occurrence: the
  first occurrence in one sequence corresponds to the first in the other, the
  second to the second, and so on.

The examples on this page use the order of six texts in two manuscripts:

```python
>>> import seqsim
>>> ms_a = ["Prologue", "Knight", "Miller", "Reeve", "Cook", "Man of Law"]
>>> ms_b = ["Prologue", "Knight", "Man of Law", "Miller", "Reeve", "Cook"]
```

## Ulam

`order.ulam_dist(x, y, *, normal=False)`

For permutations, the Ulam distance is the minimum number of items that must
be moved (taken out and reinserted elsewhere) to transform one order into the
other, which equals the number of items minus the length of their longest
common subsequence. `seqsim` generalizes it to the minimum number of moves,
insertions, and deletions of single items:
$|x| + |y| - M - \mathrm{LCS}(x, y)$, where $M$ is the number of items the
sequences share.

```python
>>> seqsim.order.ulam_dist(ms_a, ms_b)
1.0
>>> seqsim.order.ulam_dist(ms_a, ms_b + ["Squire"])
2.0
```

Properties
: A true distance (it is an edit distance with symmetric, unit-cost
  operations). With `normal=True`, it is divided by $|x| + |y| - M$, its
  maximum.

When to use
: The simplest measure of reordering, modelling a scribe or compiler moving
  single texts. A good default.

References
: Ulam (1972); Aldous and Diaconis (1999).

## Kendall tau

`order.kendall_tau_dissim(x, y, *, p=0.5, normal=False)`

For permutations, the Kendall tau distance is the number of pairs of items in
a different relative order, which is also the minimum number of swaps of
adjacent items needed to transform one order into the other. For sequences
with different items, `seqsim` follows the generalization of Fagin et al.
(2003) for "top-k lists", where an item missing from a sequence is considered
to come after all of its items; a pair of items both absent from one sequence
has an unknown relative order, and costs `p`. An end marker shared by both
sequences makes each item present in only one sequence count.

```python
>>> seqsim.order.kendall_tau_dissim(ms_a, ms_b)
3.0
```

Here "Man of Law" moved over three texts, so three pairs are in a different
order.

Properties
: A true distance for permutations. For sequences with different items it is
  a "near metric" but does not satisfy the triangle inequality, for any value
  of `p` (Fagin et al., 2003). With `normal=True`, it is divided by its
  maximum over all orders of the same contents; for permutations of $n$
  items this is $n(n-1)/2$, the classic normalized Kendall distance, so that
  two reversed orders score 1.0.

```python
>>> seqsim.order.kendall_tau_dissim(ms_a, ms_b, normal=True)
0.2
>>> seqsim.order.kendall_tau_dissim(ms_a, ms_a[::-1], normal=True)
1.0
```

When to use
: When a text moved far should count more than a text moved to the next
  position.

References
: Kendall (1938); Fagin et al. (2003).

## Kendall correlation of the shared items

`order.kendall_tau_simil(x, y, *, normal=False)`

Kendall's rank correlation $\tau$ between the orders of the items both
sequences share, ignoring all other items: $1 - 4d / (n(n-1))$, where $n$ is
the number of shared items and $d$ the number of pairs of them in a
different relative order. It is 1.0 when the shared items are in the same
order and -1.0 when they are reversed; for permutations it equals Kendall's
tau-a and tau-b, which coincide without ties. It asks only whether the
material two witnesses have in common is in the same order, independently of
how much each of them selected.

```python
>>> seqsim.order.kendall_tau_simil(ms_a, ms_b)
0.6
>>> seqsim.order.kendall_tau_simil(ms_a, ["Squire", "Man of Law", "Knight"])
-1.0
```

With fewer than two shared items the correlation is undefined, and `nan` is
returned. A minimum number of shared items, below which the correlation is
not worth reporting, is left to the caller. With `normal=True`, the
correlation is mapped to range [0..1] as $(1 + \tau) / 2$.

References
: Kendall (1938).

## Spearman footrule

`order.footrule_dissim(x, y, *, ell=None, normal=False)`

The sum, over all items, of the absolute difference of their positions in
both sequences. Items missing from a sequence are placed at position `ell`,
by default one after the end of the longest sequence. The footrule is always
within a factor of two of the Kendall tau distance (Diaconis and Graham,
1977).

```python
>>> seqsim.order.footrule_dissim(ms_a, ms_b)
6.0
```

Properties
: A true distance for permutations, and for any sequences when the same
  fixed `ell`, larger than the length of every sequence being compared, is
  used for all comparisons (it is then the L1 distance between vectors of
  positions). With the default `ell`, which depends on the pair of sequences,
  the triangle inequality can fail.

When to use
: As a simple and fast alternative to Kendall tau. Pass a fixed `ell` when
  building a distance matrix.

References
: Spearman (1906); Diaconis and Graham (1977); Fagin et al. (2003).

## Cayley

`order.cayley_dissim(x, y, *, normal=False)`

For permutations, the Cayley distance is the minimum number of exchanges of
two items (not necessarily adjacent) needed to transform one order into the
other, that is, the number of items minus the number of cycles of the
permutation relating them. Items present in only one sequence add one each.

```python
>>> seqsim.order.cayley_dissim(ms_a, ms_b)
3.0
```

Properties
: A true distance for permutations, but not for sequences with repeated
  items (e.g., `"aab"`, `"caab"`, and `"baac"`).

References
: Cayley (1849); Diaconis (1988).

## Block interchange

`order.block_interchange_dissim(x, y, *, normal=False)`

The minimum number of exchanges of two blocks of consecutive items (of any
size, not necessarily adjacent) needed to transform one order into the other,
computed with the formula of Christie (1996), $(n + 1 - c) / 2$, where $c$ is
the number of cycles of the "cycle graph" of the permutation. Moving a single
block elsewhere is a special case. Items present in only one sequence add one
each.

```python
>>> seqsim.order.block_interchange_dissim(ms_a, ms_b)
1.0
>>> quires = list("ABCDEFGH")
>>> seqsim.order.block_interchange_dissim(quires, list("EFGHABCD"))
1.0
```

Properties
: A true distance for permutations. For sequences with different or repeated
  items, no counterexample to the triangle inequality is known, but it is not
  proven.

When to use
: When groups of texts (quires, gatherings, booklets) are bound or copied in
  a different order.

References
: Christie (1996).

## Breakpoints

`order.breakpoint_dissim(x, y, *, normal=False)`

An adjacency is a pair of consecutive items, including a start and an end
marker, so that a sequence of $n$ items has $n + 1$ adjacencies. The
breakpoint dissimilarity is half the size of the symmetric difference of the
adjacencies of both sequences; for permutations, it is the number of
adjacencies of one sequence that are broken in the other. Adjacencies are
ordered, as the direction of reading matters for texts.

```python
>>> seqsim.order.breakpoint_dissim(ms_a, ms_b)
3.0
```

Moving "Man of Law" broke three adjacencies of the first manuscript:
"Knight" followed by "Miller", "Cook" followed by "Man of Law", and "Man of
Law" at the end.

With `boundaries=False`, only the $n - 1$ adjacencies between items are
counted. This is the right choice for fragments, which can begin and end
anywhere, so that their first and last texts say nothing about order:

```python
>>> seqsim.order.breakpoint_dissim(ms_a, ms_b, boundaries=False)
2.0
```

`order.breakpoint_simil(x, y, *, boundaries=True)` is the share of
adjacencies preserved, one minus the normalized dissimilarity, computed
directly so that its value is exact. Applied to the sequences reduced to
their shared items (see below), and without boundaries, it measures how well
two witnesses agree on the succession of the material they share, a finer
measure than Kendall's correlation when only a few blocks have moved:

```python
>>> shared_a, shared_b = seqsim.order.restrict_to_shared(ms_a, ms_b + ["Squire"])
>>> seqsim.order.breakpoint_simil(shared_a, shared_b, boundaries=False)
0.6
```

Properties
: Symmetric and satisfies the triangle inequality, with or without
  boundaries, but with repeated items different sequences can have the same
  adjacencies (e.g., `"abacada"` and `"acabada"`). Without boundaries,
  sequences of fewer than two items have no adjacencies and score 0.0
  against each other (e.g., `"a"` and `"b"`).

When to use
: A robust, local measure of reordering, used by Spencer et al. (2003) for
  the order of the *Canterbury Tales*. It saturates when there were many
  rearrangements; see IEBP below.

References
: Sankoff and Blanchette (1998); Spencer et al. (2003).

## IEBP

`order.iebp_estimate(x, y, *, boundaries=True, normal=False)`

The breakpoint dissimilarity underestimates the number of rearrangements
when there were many, as later rearrangements can break adjacencies that were
already broken. The IEBP method ("Inverse of the Expected BreakPoint
distance") of Wang and Warnow (2001) estimates the number of rearrangements
$k$ whose expected number of breakpoints is closest to the observed one.
`seqsim` uses the formulas of Spencer et al. (2003) for linear orders
rearranged by transpositions (moves of blocks of one or more items). Only the
items shared by both sequences are considered.

```python
>>> seqsim.order.iebp_estimate(ms_a, ms_b)
1.0
>>> seqsim.order.iebp_estimate(ms_a, ms_b, normal=True)
0.16666666666666666
```

With `normal=True`, the estimate is divided by the number of shared items, as
done by Spencer et al. (2003); the result is not bounded by one. With
`boundaries=False`, breakpoints are counted only between items, and the
expected number is the sum of the terms of the formulas for these interior
positions. Note that a block moved to the beginning or the end breaks a
single interior adjacency, fewer than the average transposition, and may be
estimated as zero moves.

Properties
: An estimator, not a measure of distance: sequences whose shared items are
  in the same order have an estimate of zero regardless of any other items,
  so it is not available through `distance()`. When the observed breakpoints
  are close to the maximum expected under the model (with few items or many
  rearrangements), the estimate saturates and should be read as "many".

References
: Wang and Warnow (2001); Spencer et al. (2003).

## Restricting to the shared items

`order.restrict_to_shared(x, y, *, repeats="occurrence")`

Reduces both sequences to the items they share, each in its own order. Any
measure can then compare the order of the shared material alone, for
example of two witnesses that selected different texts from the same
collection:

```python
>>> seqsim.order.restrict_to_shared(["a", "b", "X", "c"], ["c", "Y", "a", "b"])
(['a', 'b', 'c'], ['c', 'a', 'b'])
```

Repeated items are matched by occurrence, as in all measures of the module
(`repeats="occurrence"`); with `repeats="first"`, each sequence is first
reduced to the first occurrence of each item.

```python
>>> seqsim.order.restrict_to_shared("abab", "bab")
(['a', 'b', 'b'], ['b', 'a', 'b'])
>>> seqsim.order.restrict_to_shared("abab", "bab", repeats="first")
(['a', 'b'], ['b', 'a'])
```
