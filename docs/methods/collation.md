# Collations

The `seqsim.collation` module works on a **collation**: witnesses aligned
into variant locations (**loci**), as an alignment or a critical apparatus
gives them. A collation maps each siglum to a sequence of readings of the
same length, one per locus. A reading is any hashable value (a word, a
normalized spelling, the identifier of a group of readings), and a marker of
omission such as `""` is a reading like any other. `None` marks missing data:
a lacuna, damage, or a passage the witness never had.

Unlike the measures of the other modules, which compare two sequences on
their own, the measures of this module need the whole collation, because
they weigh what two witnesses share by how common it is among all of them.

```python
>>> import seqsim
>>> from seqsim import collation
>>> witnesses = {
...     "A": ["in", "principio", "erat", "verbum"],
...     "B": ["in", "principio", "erat", "sermo"],
...     "C": ["in", "initio", "erat", "sermo"],
...     "D": ["in", "initio", "fuit", None],
... }
```

## Bayesian affinity

`collation.affinity(collation, *, witnesses=None, weight=None,
counting="locus", prior=(1, 9), level=0.95, grid=None)`

Estimates, for every pair of witnesses, an **affinity** rho: the probability
that one witness takes its reading from the line of the other, rather than
from the tradition at large. Agreement in a rare reading is strong evidence
of affinity, agreement in a reading most witnesses share is weak evidence,
and every disagreement counts against it. The result is a posterior mean
with a credible interval, together with the raw counts:

```python
>>> result = collation.affinity(witnesses)
>>> ab = result["A", "B"]
>>> round(ab.rho, 4), round(ab.low, 4), round(ab.high, 4)
(0.1065, 0.0031, 0.3464)
>>> ab.comparable, ab.agreements, ab.agreement
(4, 3, 0.75)
```

The measure is the "Bayesian affinity ρ" of the *Proximity according to
shared readings* tool of the
[Corpus Avesticum Berolinense](https://cab.geschkult.fu-berlin.de/) (CAB,
Freie Universität Berlin; Alberto Cantera), whose method is documented only
by the source code of the tool (version 34, consulted 24 November 2026).
`seqsim` reproduces CAB's computation exactly when asked to (see
[below](#reproducing-cab)), and by default makes two changes, explained
below: it counts one observation per locus, and it integrates the posterior
accurately.

### The model

At a locus where witness A has a reading held by `k` of the `N` witnesses
with a reading there, another witness drawn at random has the same reading
with probability `q = (k - 1) / (N - 1)`. The model takes witness B to have
A's reading because it shares A's line with probability rho, and otherwise to
have a reading drawn from the tradition, which agrees with A with
probability `q`:

    P(agreement)    = rho + (1 - rho) q
    P(disagreement) = (1 - rho) (1 - q)

The factor `(1 - q)` does not depend on rho and drops out. With a Beta(a, b)
prior, the posterior density is proportional to

    rho^(a - 1) (1 - rho)^(b - 1 + D) prod_j (q_j + (1 - q_j) rho)^(w_j)

where the product runs over the shared readings, each with its weight
`w_j`, and `D` is the weighted number of disagreements. The default prior,
Beta(1, 9), is CAB's: it has mean 0.1 and expects affinities to be small.

`q` is counted among all the witnesses of the collation with a reading at
the locus, also when `witnesses` restricts the pairs that are computed.
Loci where fewer than two witnesses have a reading are ignored, and so are,
for a pair, the loci where either witness has missing data.

### Counting the evidence

With `counting="locus"` (the default), every locus where both witnesses have
a reading is one observation. The model is directional (it takes `q` from
the reading of A), and the measure is made symmetric by averaging the
log-likelihood of both directions: an agreement contributes
`w log(q + (1 - q) rho)`, and a disagreement `(w_A + w_B) / 2 log(1 - rho)`.

With `counting="reading"`, the counting of CAB, every reading held by at
least one of the two witnesses at a locus is an observation: "shared" if
both hold it, a "mismatch" otherwise. An agreement then counts once, and a
disagreement twice, once for each of the two readings:

```python
>>> loc = collation.affinity(witnesses)["A", "D"]
>>> rdg = collation.affinity(witnesses, counting="reading")["A", "D"]
>>> loc.differing, rdg.differing
(2.0, 4.0)
```

### Weights

`weight` is a function of the index of a locus and a reading, returning the
weight of that reading (by default 1). On a critical apparatus, CAB gives
explicit variants weight 1, other recorded readings 0.25 (or 0.5, depending
on the selection), and the reading of the edited text (the lemma) weight 0,
so that agreement in the edited text is no evidence of affinity, while
agreement in a variant is:

```python
>>> lemma = witnesses["A"]
>>> result = collation.affinity(
...     witnesses, weight=lambda i, reading: 0.0 if reading == lemma[i] else 1.0
... )
>>> result["A", "B"].shared, result["A", "B"].differing
(0.0, 0.5)
```

### The integration

CAB computes the posterior on a grid of 301 equally spaced values of rho.
The grid cannot resolve a posterior narrower than its spacing: with many
disagreements, as on a long apparatus, the posterior lies below the first
grid point, 1/602, and the reported mean and both ends of the interval are
that point, whatever the evidence. By default, `seqsim` integrates on the
logit scale, over a window adapted to the posterior, with an accuracy better
than one part in a million against closed forms:

```python
>>> collation.posterior([(0.01, 30.0)], 20000.0, grid=301)
(0.0016611295681063123, 0.0016611295681063123, 0.0016611295681063123)
>>> f"{collation.posterior([(0.01, 30.0)], 20000.0)[0]:.3g}"
'5.86e-05'
```

`collation.posterior(shared, differing, prior=(1, 9), level=0.95, grid=None)`
computes the posterior from the evidence directly, given as pairs `(q, w)`
and the weighted number of disagreements.

### Reproducing CAB

`counting="reading"`, `grid=301`, and the weights of CAB reproduce its
values to the precision of floating-point arithmetic; the test suite checks
this against values computed with CAB's own functions. CAB also groups
readings after normalizing spelling and excludes some kinds of evidence
(gaps, *deest*, abbreviations); those choices belong to the preparation of
the collation, before `affinity`. CAB's "coverage-adjusted ranking", the
affinity multiplied by `sqrt(n / (n + tau))` for `n` comparable loci, is a
heuristic for ordering links in its display, not part of the posterior, and
is not offered.

### Properties and use

The affinity is a similarity in range [0..1], symmetric, and not defined for
a witness with itself; `collation.affinity_matrix` returns it as a matrix
with 1.0 on the diagonal (and NaN for pairs without evidence). It is not a
metric, and `1 - rho` is a dissimilarity that can be given to tree builders.
It is a measure of shared innovation in the spirit of the common-error
method, weighted by rarity, but it does not distinguish innovations from
retained readings: when the lemma is the archetype's reading, weighting it
0 comes close.
