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

The measure is the "Bayesian affinity ρ" of the
[*Proximity according to shared readings*](https://ada.geschkult.fu-berlin.de/manuscript_analysis/proximity.html)
tool of the [Corpus Avesticum Berolinense](https://cab.geschkult.fu-berlin.de/)
(CAB, Freie Universität Berlin), which describes it as a "provisional
stemmatic affinity measure". Its method is documented only by the source code
of the tool (`common.js` version 30, loaded by the page of version 34,
consulted 24 November 2026). By default, `seqsim` counts the evidence as CAB
does on aligned transliterations and integrates the posterior more
accurately than CAB's grid; it gives the values computed by CAB's
functions when asked to, under the conditions given
[below](#reproducing-cab).

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

CAB counts the evidence in two ways, one for each of its sources, and
`seqsim` offers both.

With `counting="locus"` (the default), as CAB does on aligned
transliterations, every locus where both witnesses have a reading is one
observation. The model is directional (it takes `q` from the reading of A),
and averaging the log-likelihood of both directions makes it symmetric: an
agreement contributes `w log(q + (1 - q) rho)`, and a disagreement
`(w_A + w_B) / 2 log(1 - rho)` (CAB gives every reading weight 1).

With `counting="reading"`, as CAB does on a critical apparatus, every
reading held by at least one of the two witnesses at a locus is an
observation: "shared" if both hold it, a "mismatch" otherwise. An agreement
then counts once, and a disagreement twice, once for each of the two
readings, so that the affinities are lower:

```python
>>> loc = collation.affinity(witnesses)["A", "D"]
>>> rdg = collation.affinity(witnesses, counting="reading")["A", "D"]
>>> loc.differing, rdg.differing
(2.0, 4.0)
```

### Weights

`weight` is a function of the index of a locus and a reading, returning the
weight of that reading (by default 1). On a critical apparatus, CAB gives
explicit variants and omissions (types "var", "´var", "om", "om.") weight
1, readings of type "reading" 0.25 (0 when only variants are selected),
readings of other types 0.5 (only when all evidence is selected), and the
reading of the edited text (the lemma) weight 0; a group of readings merged
from several types takes the largest of their weights. Agreement in the
edited text is then no evidence of affinity, while agreement in a variant
is:

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
For most collations the grid changes the mean little, and an interval whose
lower end is the first grid point, 1/602, usually has a true lower end just
below it. The grid cannot resolve a posterior narrower than its spacing,
though: when the whole posterior lies below 1/602, as it can on very long
collations where every reading counts, the reported mean and both ends of
the interval are that point, whatever the evidence, and pairs tie. By
default, `seqsim` integrates on the logit scale, over a window adapted to
the posterior, with the tails beyond it integrated exactly. Against
high-precision quadrature the mean agrees to a relative error of about
1e-10 and the ends of the interval to about 2e-5; the tests check the mean
against closed forms to 1e-7 and the interval to 1e-5:

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

With `grid=301`, `affinity` caps `q` at 0.999999, as CAB's counting does,
and gives the values of ρ and of its interval computed by CAB's functions
to the precision of floating-point arithmetic; the test suite checks both
modes against values computed with CAB's own functions on synthetic
collations of eight witnesses. No output of CAB's tool on real data has
been compared. (`posterior(..., grid=301)` reproduces CAB's
`bayesPosterior`, which does not cap `q`.)

- **Transliteration mode** (`pairMsStats`): the default counting with
  `grid=301`, on a collation holding only the witnesses selected in CAB,
  since CAB counts `q` among them.
- **Apparatus mode** (`appPairMatrix`): `counting="reading"`, `grid=301`, and
  CAB's weights, with `witnesses` for the selection (`q` is counted among
  all witnesses), under three conditions. Every witness must have one
  reading per locus: CAB merges corrections such as *a.c.* and *p.c.* into
  one siglum, which can then hold two readings, and `seqsim` cannot
  represent that. The weight of a reading must be that of CAB's group of
  readings, the largest over the types of the readings merged into it. The
  readings CAB discards (gaps, *deest*, abbreviations, *bis* witnesses) must
  be removed beforehand, as missing data.

Two limits of CAB's implementation are not reproduced. In transliteration
mode, CAB stores the counts of readings in 8-bit integers, which wrap
around with 256 or more witnesses. It also groups shared readings by `q`
rounded to eight decimals, which only matters if two different values of
`q` agree to that precision. The fixtures are collations of eight
witnesses. In apparatus mode, CAB also counts the comparable loci and the
agreements of its raw agreement only at loci with a variant of the selected
kind, while `comparable` and `agreements` count every locus where both
witnesses have a reading.

CAB's normalization of spelling belongs to the preparation of the
collation, before `affinity`. CAB's "coverage-adjusted ranking", the
affinity multiplied by `sqrt(n / (n + tau))` for `n` comparable loci, is a
heuristic for ordering links in its display, not part of the posterior, and
is not offered.

### Properties and use

The affinity is a similarity in range [0..1], symmetric, and not defined for
a witness with itself; `collation.affinity_matrix` returns it as a matrix
with 1.0 on the diagonal (and NaN for pairs without evidence). It is not a
metric. `1 - rho` is a dissimilarity, but it saturates for distant pairs,
whose affinities are all close to 0; `-log rho` would be additive along the
paths of a tree if affinities multiplied along its branches, which is a
motivation for using it with tree builders such as neighbor-joining, not
a property the estimated affinities are known to have.
It is a rarity-weighted measure of agreement. It does not distinguish
innovations from retained readings, nor the direction of copying, and
giving the lemma weight 0 does not make the remaining agreements shared
innovations.

### What the interval means

The likelihood treats the loci as independent observations and `q` as
known, takes `q` from all witnesses including the pair itself (so that a
reading held only by the two witnesses still has `q > 0`), and uses weights
as exponents. The loci of a tradition are not independent and its witnesses
are not exchangeable draws, so the credible interval summarizes this
working likelihood; it is not a calibrated statement of uncertainty about
the relation of two witnesses. The prior matters: CAB's Beta(1, 9) has
mean 0.1, above most affinities in long collations, and with only two
witnesses at a locus every agreement has `q = 1` and is no evidence, while
a disagreement always counts. The averaging of the two directions in
`counting="locus"` is a device of `seqsim` for unequal weights; with unit
weights, as in CAB's transliteration mode, it changes nothing.
