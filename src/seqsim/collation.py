"""
Measures over a collation: witnesses aligned into variant locations.

A collation maps each witness (a siglum) to a sequence of readings of the
same length, one per variant location ("locus"), as an alignment or a
critical apparatus gives them. A reading is any hashable value: a word, a
normalized spelling, an identifier of a reading group, or a marker of
omission (such as `""`), which counts as a reading. `None` marks missing
data (a lacuna, a damaged or unread passage): the witness says nothing about
that locus. Missing data must be `None`: any other value, including a float
NaN, is a reading.

Unlike the measures of the other modules, which compare two sequences on
their own, the measures of this module need the whole collation, because
they weigh the readings two witnesses share by how common those readings
are among all the witnesses.

Bayesian affinity
-----------------

`affinity` estimates, for every pair of witnesses, a "Bayesian affinity"
rho: the probability that one witness takes a reading from the line of the
other rather than from the tradition at large. It follows the "Bayesian
affinity" of the *Proximity according to shared readings* tool of the Corpus
Avesticum Berolinense (CAB, Freie Universität Berlin), which calls it a
"provisional stemmatic affinity measure". Its definition is given only by the
source code of that tool (`common.js` version 30, loaded by the page of
version 34, at https://ada.geschkult.fu-berlin.de/manuscript_analysis/,
consulted 24 November 2026).

At a locus where witness A has a reading shared by `k` of the `N` witnesses
with a reading there, a witness B drawn at random from the others has the
same reading with probability `q = (k - 1) / (N - 1)`. The model mixes two
sources: with probability rho, B has A's reading because it shares A's line;
otherwise it has a reading drawn from the tradition, which agrees with A
with probability `q`. Hence

    P(agreement)    = rho + (1 - rho) q = q + (1 - q) rho
    P(disagreement) = (1 - rho) (1 - q)

An agreement in a rare reading (small `q`) is strong evidence of affinity,
an agreement in a reading most witnesses have (large `q`) is weak evidence,
and a disagreement always counts against it. The factor `(1 - q)` does not
depend on rho and is left out. With a Beta(a, b) prior (Beta(1, 9) by
default, as in CAB, which expects affinities to be small), the posterior
density is proportional to

    rho^(a - 1) (1 - rho)^(b - 1 + D) prod_j (q_j + (1 - q_j) rho)^(w_j)

where `D` is the weighted number of disagreements and the product runs over
the shared readings, each with its weight `w_j`. `affinity` reports the
posterior mean and an equal-tailed credible interval.

CAB counts the evidence in two ways, one for each of its sources, and both
are offered:

  * `counting="locus"` (the default), as CAB does on aligned
    transliterations: every locus where both witnesses have a reading is one
    observation. An agreement contributes `w log(q + (1 - q) rho)` and a
    disagreement `(w_A + w_B) / 2 log(1 - rho)`, with `w_A` and `w_B` the
    weights of the two readings (CAB gives every reading weight 1). The
    model is directional (it takes `q` from the reading of A), and averaging
    the log-likelihoods of both directions gives these terms, since two
    witnesses that agree have the same reading.
  * `counting="reading"`, as CAB does on a critical apparatus: every reading
    at the locus held by at least one of the two witnesses is an
    observation, "shared" if both hold it and a "mismatch" otherwise. An
    agreement counts once and a disagreement twice (once for each of the two
    readings).

Weights let a reading count for less or not at all. On a critical
apparatus, CAB gives explicit variants and omissions weight 1 (types "var",
"´var", "om" and "om."),
readings of type "reading" 0.25 (0 when only variants are selected),
readings of any other type 0.5 (only when all evidence is selected), and
the reading of the edited text (the lemma) weight 0, so that agreement in
the edited text is not evidence while agreement in a variant is.

CAB integrates the posterior on a grid of 301 points of rho, and caps `q`
at 0.999999. A grid cannot resolve a posterior narrower than its spacing:
when the posterior lies below the first grid point (1/602), as it can on
long collations where every reading counts, the reported mean and interval
are that point. By default, `affinity` integrates on the logit scale over a
window adapted to the posterior; `grid=301` follows the computation of
CAB, with its cap on `q`.

Reproducing CAB
***************

`grid=301` with the default counting gives the values of ρ and of its
interval computed by CAB's transliteration mode
(`pairMsStats`) when the collation holds only the witnesses selected in CAB
(CAB counts `q` among them). `grid=301`, `counting="reading"` and weights
equal to CAB's give those of its apparatus mode (`appPairMatrix`), where `q` is
counted among all witnesses and `witnesses` gives the selection, under
three conditions: every witness has one reading per locus (CAB merges
corrections such as *a.c.* and *p.c.* into one siglum, which can then hold
two readings); the weight of a reading is that of CAB's group of readings,
the largest over the types of the readings merged into it; and the readings
CAB discards (gaps, *deest*, abbreviations, *bis* witnesses) are removed
beforehand, as missing data. CAB's normalization of spelling belongs to the
preparation of the collation. Two limits of CAB's implementation are not
reproduced: in transliteration mode, CAB stores the counts of readings in
8-bit integers, which wrap around with 256 or more witnesses, and it groups
shared readings by `q` rounded to eight decimals, which only matters if two
different values of `q` agree to that precision. In apparatus mode, CAB
counts the comparable loci and the agreements of its raw agreement only at
loci with a variant of the selected kind, while `Affinity.comparable` and
`Affinity.agreements` count every locus where both witnesses have a
reading. The tests check both modes against values computed with CAB's own
functions on synthetic collations of eight witnesses; no output of CAB's
tool on real data has been compared.

What the interval means
***********************

The likelihood treats the loci as independent observations and `q` as
known, takes `q` from all witnesses including the pair itself (so that a
reading held only by the two witnesses still has `q > 0`), and uses weights
as exponents. Loci of a tradition are not independent, and witnesses are
not exchangeable draws from it, so the credible interval is a summary of
this working likelihood, not a calibrated statement of uncertainty about
the relation of two witnesses. The averaging of the two directions in
`counting="locus"` is a device of `seqsim` for unequal weights; with unit
weights, as in CAB's transliteration mode, it changes nothing.

Example
*******

.. code-block:: python

    >>> collation = {
    ...     "A": ["in", "principio", "erat", "verbum"],
    ...     "B": ["in", "principio", "erat", "sermo"],
    ...     "C": ["in", "initio", "erat", "sermo"],
    ...     "D": ["in", "initio", "fuit", None],
    ... }
    >>> result = seqsim.collation.affinity(collation)
    >>> round(result["A", "B"].rho, 4), round(result["B", "C"].rho, 4)
    (0.1065, 0.1009)

A and B agree in "principio", a reading of two of the four witnesses; B
and C agree in "sermo", a reading of two of the three witnesses that
preserve the last word, which is weaker evidence.
"""

# Import Python standard libraries
import math
from collections import Counter
from collections.abc import Callable, Hashable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field

# Missing data in a collation
MISSING = None

# Highest value of `q` in the counting of CAB
_CAB_MAX_Q = 0.999999

# Points of the initial scan and of each refinement, on the logit scale
_SCAN = 241
_FINE = 400  # intervals (even, for Simpson's rule)
_MAX_STEP = 0.02
_LOGIT_RANGE = 30.0
# Where the window stops on the logit scale and the tails are integrated
# analytically: beyond these points the log density is linear in logit(rho)
_TAIL_UPPER = 40.0
_TAIL_LOWER = -40.0
_TAIL_MARGIN = 30.0
# Log-density below the maximum at which the window is cut
_CUT = 46.0

__all__ = ["MISSING", "Affinity", "affinity", "affinity_matrix", "posterior"]


@dataclass(frozen=True)
class Affinity:
    """
    The Bayesian affinity of a pair of witnesses (see `affinity`).

    :ivar rho: The posterior mean of the affinity.
    :ivar low: The lower end of the credible interval.
    :ivar high: The upper end of the credible interval.
    :ivar shared: The weighted number of shared readings.
    :ivar differing: The weighted number of disagreements.
    :ivar comparable: The number of loci where both witnesses have a reading.
    :ivar agreements: The number of those loci where they have the same
        reading.
    """

    rho: float
    low: float
    high: float
    shared: float
    differing: float
    comparable: int
    agreements: int

    @property
    def agreement(self) -> float | None:
        """The share of comparable loci with the same reading (raw agreement)."""
        return self.agreements / self.comparable if self.comparable else None


@dataclass
class _Evidence:
    """The evidence about a pair of witnesses, accumulated over the loci."""

    comparable: int = 0
    agreements: int = 0
    terms: dict[float, float] = field(default_factory=dict)  # weight by q
    shared: float = 0.0
    differing: float = 0.0


def _check_prior(prior: tuple[float, float]) -> tuple[float, float]:
    a, b = prior
    if not (a > 0 and b > 0 and math.isfinite(a) and math.isfinite(b)):
        raise ValueError("The parameters of the prior must be positive and finite.")
    return float(a), float(b)


def _log_density(
    rho: float, a: float, b: float, terms: list[tuple[float, float]], differing: float
) -> float:
    """The unnormalized log posterior density of rho (on the rho scale)."""
    value = (a - 1.0) * math.log(rho) + (b - 1.0 + differing) * math.log1p(-rho)
    for q, w in terms:
        value += w * math.log(q + (1.0 - q) * rho)
    return value


def _sigmoid(u: float) -> float:
    if u >= 0:
        return 1.0 / (1.0 + math.exp(-u))
    e = math.exp(u)
    return e / (1.0 + e)


def _log_density_logit(
    u: float, a: float, b: float, terms: list[tuple[float, float]], differing: float
) -> float:
    """The unnormalized log posterior density of logit(rho)."""
    rho = _sigmoid(u)
    # log rho and log(1 - rho), stable for large |u|
    # (each computed directly: log_rho - u would cancel when u is very
    # negative, which matters for many disagreements)
    if u >= 0:
        log_rho = -math.log1p(math.exp(-u))
        log_one_minus = -u - math.log1p(math.exp(-u))
    else:
        log_rho = u - math.log1p(math.exp(u))
        log_one_minus = -math.log1p(math.exp(u))
    value = a * log_rho + (b + differing) * log_one_minus
    for q, w in terms:
        # q = 0 gives w log(rho), which must not go through rho (it underflows)
        value += w * (log_rho if q == 0.0 else math.log(q + (1.0 - q) * rho))
    return value


def _intervals(lo: float, hi: float) -> int:
    """An even number of intervals, at least `_FINE` and no wider than `_MAX_STEP`."""
    intervals = max(_FINE, math.ceil((hi - lo) / _MAX_STEP))
    return intervals + intervals % 2


def _window(
    logd: Callable[[float], float],
    points: list[float],
    values: list[float],
    step: float,
    limits: tuple[float, float],
) -> tuple[float, float]:
    """
    The interval around the maximum where the log density is within `_CUT`.

    The interval is extended beyond the points when the log density has not
    fallen by `_CUT` at their ends (heavy tails, as with a prior `a < 1`),
    but not beyond `limits`, past which the tails are integrated
    analytically.
    """
    top = max(range(len(values)), key=values.__getitem__)
    peak = values[top]
    lo = top
    while lo > 0 and values[lo] > peak - _CUT:
        lo -= 1
    hi = top
    while hi < len(values) - 1 and values[hi] > peak - _CUT:
        hi += 1
    # One more step on each side, in case the peak falls between points
    first, last = limits
    low = max(points[lo] - step, first)
    high = min(points[hi] + step, last)
    width = max(step, 1.0)
    while low > first and logd(low) > peak - _CUT:
        low = max(low - width, first)
        width *= 2.0
    width = max(step, 1.0)
    while high < last and logd(high) > peak - _CUT:
        high = min(high + width, last)
        width *= 2.0
    return low, high


def _summary(
    points: list[float],
    density: list[float],
    values: list[float],
    level: float,
    slopes: tuple[float, float],
) -> tuple[float, float, float]:
    """
    Mean and equal-tailed interval from a density on an even grid of logit(rho).

    Integrals use Simpson's rule; the cumulative integral is also computed at
    the midpoint of every pair of intervals, and quantiles are found by
    inverting, within an interval, the integral of a linear density. Beyond
    the grid, the log density is taken as linear in logit(rho), with the
    given slopes (the rates of decay of the lower and upper tails), and the
    tails are integrated exactly.
    """
    h = points[1] - points[0]
    rate_low, rate_high = slopes
    f_low, f_high = density[0], density[-1]
    # Lower tail: density f_low e^(rate (u - u0)); rho ~ e^u there
    mass_low = f_low / rate_low
    mean_low = f_low * values[0] / (rate_low + 1.0)
    # Upper tail: density f_high e^(-rate (u - u1)); 1 - rho ~ e^(-u) there
    mass_high = f_high / rate_high
    mean_high = mass_high - f_high * (1.0 - values[-1]) / (rate_high + 1.0)
    cumulative = [mass_low]
    mean = mean_low
    for i in range(0, len(points) - 2, 2):
        f0, f1, f2 = density[i], density[i + 1], density[i + 2]
        base = cumulative[-1]
        cumulative.append(base + h / 12.0 * (5.0 * f0 + 8.0 * f1 - f2))
        cumulative.append(base + h / 3.0 * (f0 + 4.0 * f1 + f2))
        mean += (
            h / 3.0 * (f0 * values[i] + 4.0 * f1 * values[i + 1] + f2 * values[i + 2])
        )
    total = cumulative[-1] + mass_high
    mean = (mean + mean_high) / total

    def quantile(p: float) -> float:
        target = p * total
        if target < mass_low:
            return points[0] + math.log(target / mass_low) / rate_low
        if target > cumulative[-1]:
            rest = (total - target) / mass_high
            return points[-1] - math.log(max(rest, 1e-300)) / rate_high
        for j in range(1, len(cumulative)):
            if cumulative[j] >= target:
                need = target - cumulative[j - 1]
                f0, f1 = density[j - 1], density[j]
                # Solve f0 t + (f1 - f0) t^2 / (2 h) = need for t in [0..h]
                slope = (f1 - f0) / h
                if abs(slope) * h < 1e-12 * max(f0, f1, 1e-300):
                    t = need / f0 if f0 > 0 else 0.0
                else:
                    disc = max(f0 * f0 + 2.0 * slope * need, 0.0)
                    t = (math.sqrt(disc) - f0) / slope
                return points[j - 1] + min(max(t, 0.0), h)
        return points[-1]

    tail = (1.0 - level) / 2.0
    return mean, quantile(tail), quantile(1.0 - tail)


def posterior(
    shared: Iterable[tuple[float, float]],
    differing: float,
    prior: tuple[float, float] = (1.0, 9.0),
    level: float = 0.95,
    grid: int | None = None,
) -> tuple[float, float, float]:
    """
    Returns the posterior mean and credible interval of a Bayesian affinity.

    The posterior density of rho is proportional to
    `rho^(a-1) (1-rho)^(b-1+differing) prod (q + (1-q) rho)^w` over the
    shared readings `(q, w)` (see the module documentation).

    Example
    ********

    .. code-block:: python

        >>> mean, low, high = seqsim.collation.posterior([(0.1, 3.0)], 2.0)
        >>> round(mean, 4), round(low, 4), round(high, 4)
        (0.1949, 0.0166, 0.4562)

    :param shared: The shared readings, as pairs of `q` (the chance of
        agreeing in the reading by drawing from the tradition, in range
        [0..1]) and a non-negative weight.
    :param differing: The weighted number of disagreements.
    :param prior: The parameters `(a, b)` of the Beta prior. Defaults to
        `(1, 9)`.
    :param level: The probability of the credible interval. Defaults to
        0.95.
    :param grid: If given, integrate as CAB does, on this number of
        equally spaced points of rho (CAB uses 301), and take the ends of the
        interval as the first grid points where the cumulative probability
        reaches each tail. Defaults to `None`, an integration on the logit
        scale adapted to the posterior.
    :return: The posterior mean, the lower end, and the upper end of the
        credible interval.
    """

    a, b = _check_prior(prior)
    if not 0.0 < level < 1.0:
        raise ValueError("The level of the interval must be in range (0..1).")
    if not (math.isfinite(differing) and differing >= 0):
        raise ValueError("The number of disagreements must be finite and non-negative.")
    terms = []
    for q, w in shared:
        if not 0.0 <= q <= 1.0:
            raise ValueError("Every q must be in range [0..1].")
        if not (math.isfinite(w) and w >= 0):
            raise ValueError("Weights must be finite and non-negative.")
        if w:
            terms.append((float(q), float(w)))

    if grid is not None:
        return _posterior_grid(a, b, terms, float(differing), level, grid)

    def logd(u: float) -> float:
        return _log_density_logit(u, a, b, terms, float(differing))

    # Beyond these limits the log density is linear in u = logit(rho), with
    # slope a + (weights of q = 0) below and -(b + differing) above, as long
    # as rho is negligible against every positive q
    positive = [q for q, _ in terms if q > 0.0]
    # (the posterior lies around rho ~ 1 / (b + differing), so the lower
    # limit also moves with the number of disagreements)
    lower = min(_TAIL_LOWER, -math.log1p(b + float(differing)) - _TAIL_MARGIN)
    if positive:
        lower = min(lower, math.log(min(positive)) - _TAIL_MARGIN)
    limits = (lower, _TAIL_UPPER)
    slopes = (a + sum(w for q, w in terms if q == 0.0), b + float(differing))

    # A coarse scan of the logit scale, then two refinements of the window
    step = 2.0 * _LOGIT_RANGE / (_SCAN - 1)
    points = [-_LOGIT_RANGE + i * step for i in range(_SCAN)]
    values = [logd(u) for u in points]
    lo, hi = _window(logd, points, values, step, limits)
    for _ in range(2):
        intervals = _intervals(lo, hi)
        step = (hi - lo) / intervals
        points = [lo + i * step for i in range(intervals + 1)]
        values = [logd(u) for u in points]
        new_lo, new_hi = _window(logd, points, values, step, limits)
        lo, hi = max(lo, new_lo), min(hi, new_hi)
    intervals = _intervals(lo, hi)
    step = (hi - lo) / intervals
    points = [lo + i * step for i in range(intervals + 1)]
    values = [logd(u) for u in points]
    peak = max(values)
    density = [math.exp(v - peak) for v in values]
    rhos = [_sigmoid(u) for u in points]
    mean, u_low, u_high = _summary(points, density, rhos, level, slopes)
    return mean, _sigmoid(u_low), _sigmoid(u_high)


def _posterior_grid(
    a: float,
    b: float,
    terms: list[tuple[float, float]],
    differing: float,
    level: float,
    steps: int,
) -> tuple[float, float, float]:
    """The integration of CAB: a midpoint grid on the rho scale."""
    if steps < 2:
        raise ValueError("The grid needs at least two points.")
    rhos = [(i + 0.5) / steps for i in range(steps)]
    logs = [_log_density(r, a, b, terms, differing) for r in rhos]
    peak = max(logs)
    weights = [math.exp(v - peak) for v in logs]
    total = sum(weights)
    mean = sum(w * r for w, r in zip(weights, rhos)) / total
    tail = (1.0 - level) / 2.0
    cumulative, low, high, low_set = 0.0, rhos[0], rhos[-1], False
    for w, r in zip(weights, rhos):
        cumulative += w / total
        if not low_set and cumulative >= tail:
            low, low_set = r, True
        if cumulative >= 1.0 - tail:
            high = r
            break
    return mean, low, high


def _loci(
    collation: Mapping[Hashable, Sequence[Hashable]],
) -> tuple[list[Hashable], int]:
    sigla = list(collation)
    if len(sigla) < 2:
        raise ValueError("A collation needs at least two witnesses.")
    lengths = {len(collation[s]) for s in sigla}
    if len(lengths) != 1:
        raise ValueError("All witnesses of a collation must have the same length.")
    return sigla, lengths.pop()


def affinity(
    collation: Mapping[Hashable, Sequence[Hashable]],
    *,
    witnesses: Iterable[Hashable] | None = None,
    weight: Callable[[int, Hashable], float] | None = None,
    counting: str = "locus",
    prior: tuple[float, float] = (1.0, 9.0),
    level: float = 0.95,
    grid: int | None = None,
) -> dict[tuple[Hashable, Hashable], Affinity]:
    """
    Returns the Bayesian affinity of every pair of witnesses of a collation.

    See the module documentation for the model. The frequency of a reading
    at a locus is counted among all the witnesses of the collation with a
    reading there, also when `witnesses` restricts the pairs computed.
    Loci where fewer than two witnesses have a reading are ignored.

    Example
    ********

    .. code-block:: python

        >>> collation = {
        ...     "A": ["a", "b", "c", "x"],
        ...     "B": ["a", "b", "c", "x"],
        ...     "C": ["a", "y", "z", "w"],
        ... }
        >>> result = seqsim.collation.affinity(collation)
        >>> result["A", "B"].agreements, result["A", "C"].agreements
        (4, 1)
        >>> result["A", "B"].rho > result["A", "C"].rho
        True

    :param collation: A mapping from sigla to sequences of readings of the
        same length; `None` marks missing data.
    :param witnesses: The sigla whose pairs are computed. Defaults to all
        witnesses.
    :param weight: A function of the index of the locus and a reading
        returning the weight of that reading (for example 0 for the reading
        of the edited text). Defaults to weight 1 for every reading.
    :param counting: `"locus"` (one observation per locus, the default, as
        CAB counts on aligned transliterations) or `"reading"` (one per
        reading of either witness, as CAB counts on a critical apparatus).
    :param prior: The parameters `(a, b)` of the Beta prior. Defaults to
        `(1, 9)`, as in CAB.
    :param level: The probability of the credible interval. Defaults to
        0.95.
    :param grid: The number of points of CAB's grid integration (301 in
        CAB), which also caps `q` at 0.999999 as CAB does, or `None` (the
        default) for an adaptive integration; see `posterior`.
    :return: A dictionary from pairs of sigla, in both orders, to their
        `Affinity`.
    """

    if counting not in ("locus", "reading"):
        raise ValueError('`counting` must be "locus" or "reading".')
    _check_prior(prior)
    if not 0.0 < level < 1.0:
        raise ValueError("The level of the interval must be in range (0..1).")
    if grid is not None and grid < 2:
        raise ValueError("The grid needs at least two points.")
    sigla, length = _loci(collation)
    chosen = sigla if witnesses is None else list(witnesses)
    if len(set(chosen)) != len(chosen):
        raise ValueError("`witnesses` cannot repeat a witness.")
    unknown = [s for s in chosen if s not in collation]
    if unknown:
        raise ValueError(f"Unknown witnesses: {unknown!r}.")
    pairs = [(x, y) for i, x in enumerate(chosen) for y in chosen[i + 1 :]]

    records = {pair: _Evidence() for pair in pairs}
    reading_weight = (lambda i, r: 1.0) if weight is None else weight

    for i in range(length):
        readings = {s: collation[s][i] for s in sigla}
        counts = Counter(r for r in readings.values() if r is not MISSING)
        n = sum(counts.values())
        if n < 2:
            continue
        weights = {r: float(reading_weight(i, r)) for r in counts}
        if not all(math.isfinite(w) and w >= 0 for w in weights.values()):
            raise ValueError("Weights must be finite and non-negative.")
        for pair in pairs:
            ra, rb = readings[pair[0]], readings[pair[1]]
            if ra is MISSING or rb is MISSING:
                continue
            rec = records[pair]
            rec.comparable += 1
            if ra == rb:
                rec.agreements += 1
                w = weights[ra]
                if w:
                    q = (counts[ra] - 1) / (n - 1)
                    if grid is not None:
                        q = min(q, _CAB_MAX_Q)
                    rec.terms[q] = rec.terms.get(q, 0.0) + w
                    rec.shared += w
            elif counting == "locus":
                rec.differing += (weights[ra] + weights[rb]) / 2.0
            else:
                rec.differing += weights[ra] + weights[rb]

    result = {}
    for (x, y), rec in records.items():
        if rec.comparable == 0 or (rec.shared == 0 and rec.differing == 0):
            value = Affinity(
                math.nan, math.nan, math.nan, 0.0, 0.0, rec.comparable, rec.agreements
            )
        else:
            mean, low, high = posterior(
                rec.terms.items(), rec.differing, prior, level, grid
            )
            value = Affinity(
                mean,
                low,
                high,
                rec.shared,
                rec.differing,
                rec.comparable,
                rec.agreements,
            )
        result[x, y] = result[y, x] = value
    return result


def affinity_matrix(
    collation: Mapping[Hashable, Sequence[Hashable]], **kwargs
) -> tuple[list[Hashable], list[list[float]]]:
    """
    Returns the posterior mean affinities of a collation as a matrix.

    A convenience wrapper of `affinity`, accepting the same keyword
    arguments. The diagonal is 1.0, and pairs without evidence are NaN.

    :param collation: A mapping from sigla to sequences of readings.
    :return: The sigla, in the order of the collation (or of `witnesses`),
        and the square matrix of affinities.
    """

    if kwargs.get("witnesses") is not None:
        # An iterator would be consumed by `affinity`
        kwargs["witnesses"] = list(kwargs["witnesses"])
    values = affinity(collation, **kwargs)
    witnesses = kwargs.get("witnesses")
    sigla = list(collation) if witnesses is None else witnesses
    return sigla, [[1.0 if x == y else values[x, y].rho for y in sigla] for x in sigla]
