"""
test_collation
==============

Tests for the `collation` module of the `seqsim` package: the Bayesian
affinity of the *Proximity according to shared readings* tool of the Corpus
Avesticum Berolinense (CAB).

The posterior is checked against closed forms (when every shared reading has
`q = 0` or `q = 1`, it is a Beta distribution). The counting and the grid
integration of CAB are checked against values computed with the functions
`bayesPosterior` and `appPairMatrix` of CAB (`common.js`, version 34),
recorded in `tests/data/cab_affinity.json`.
"""

# Import Python standard libraries
import json
import math
import pathlib

import pytest

# Import the library being tested
from seqsim import collation
from seqsim.collation import Affinity, affinity, affinity_matrix, posterior

DATA = pathlib.Path(__file__).parent / "data" / "cab_affinity.json"


def _beta_quantile_1n(p, n):
    """Quantile of Beta(1, n)."""
    return 1.0 - (1.0 - p) ** (1.0 / n)


# The posterior


@pytest.mark.parametrize("prior", [(1, 9), (1, 1), (0.5, 0.5), (2, 3), (0.1, 0.1)])
@pytest.mark.parametrize("w", [0, 1, 5, 300, 1e4])
@pytest.mark.parametrize("d", [0, 1, 10, 1e3, 1e5, 1e7])
def test_posterior_mean_is_beta_mean(prior, w, d):
    # With q = 0, every shared reading multiplies the density by rho^w
    a, b = prior
    mean, low, high = posterior([(0.0, w)] if w else [], d, prior)
    expected = (a + w) / (a + b + w + d)
    assert mean == pytest.approx(expected, rel=1e-7)
    assert low < mean < high


@pytest.mark.parametrize("n", [1, 9, 100, 1e4, 1e6, 1e9])
def test_posterior_interval_is_beta_interval(n):
    _, low, high = posterior([], n - 1, (1, 1))
    assert low == pytest.approx(_beta_quantile_1n(0.025, n), rel=1e-5)
    assert high == pytest.approx(_beta_quantile_1n(0.975, n), rel=1e-5)
    _, low, high = posterior([], n - 1, (1, 1), level=0.5)
    assert low == pytest.approx(_beta_quantile_1n(0.25, n), rel=1e-5)
    assert high == pytest.approx(_beta_quantile_1n(0.75, n), rel=1e-5)


def test_posterior_ignores_uninformative_readings():
    # q = 1 (a reading of all witnesses) and weight 0 carry no information
    reference = posterior([], 3.0)
    assert posterior([(1.0, 5.0)], 3.0) == pytest.approx(reference, rel=1e-9)
    assert posterior([(0.2, 0.0)], 3.0) == pytest.approx(reference, rel=1e-9)


def _brute_mean(terms, d, a=1.0, b=9.0, n=200_000):
    total = weighted = 0.0
    logs = []
    for i in range(n):
        rho = (i + 0.5) / n
        logs.append(
            (a - 1) * math.log(rho)
            + (b - 1 + d) * math.log1p(-rho)
            + sum(w * math.log(q + (1 - q) * rho) for q, w in terms)
        )
    peak = max(logs)
    for i, value in enumerate(logs):
        density = math.exp(value - peak)
        total += density
        weighted += density * (i + 0.5) / n
    return weighted / total


@pytest.mark.parametrize(
    "terms,d",
    [
        ([(0.1, 3.0), (0.4, 7.5)], 12.0),
        ([(0.02, 40.0), (0.3, 10.0), (0.6, 100.0)], 150.0),
        ([(0.5, 2.0)], 0.0),
    ],
)
def test_posterior_mean_general(terms, d):
    assert posterior(terms, d)[0] == pytest.approx(_brute_mean(terms, d), rel=1e-6)


def test_posterior_more_rare_agreement_more_affinity():
    rare = posterior([(0.05, 10.0)], 20.0)[0]
    common = posterior([(0.6, 10.0)], 20.0)[0]
    more_differing = posterior([(0.05, 10.0)], 40.0)[0]
    assert rare > common
    assert rare > more_differing


@pytest.mark.parametrize(
    "kwargs",
    [
        {"prior": (0, 9)},
        {"prior": (1, -1)},
        {"prior": (1, math.inf)},
        {"level": 0},
        {"level": 1},
        {"differing": -1},
        {"shared": [(1.5, 1.0)]},
        {"shared": [(0.5, -1.0)]},
        {"grid": 1},
    ],
)
def test_posterior_errors(kwargs):
    args = {"shared": [], "differing": 1.0, **kwargs}
    with pytest.raises(ValueError):
        posterior(**args)


def test_grid_cannot_resolve_narrow_posteriors():
    # With many disagreements, the posterior lies below the first point of
    # CAB's grid (1/602): the grid returns that point for the mean and both
    # ends of the interval, whatever the evidence
    assert posterior([(0.01, 30.0)], 20000.0, grid=301) == (1 / 602, 1 / 602, 1 / 602)
    assert posterior([(0.01, 60.0)], 20000.0, grid=301) == (1 / 602, 1 / 602, 1 / 602)
    fine = posterior([(0.01, 30.0)], 20000.0)[0], posterior([(0.01, 60.0)], 20000.0)[0]
    assert fine[0] < fine[1] < 1 / 602


# The counting


COLLATION = {
    "A": ["in", "principio", "erat", "verbum"],
    "B": ["in", "principio", "erat", "sermo"],
    "C": ["in", "initio", "erat", "sermo"],
    "D": ["in", "initio", "fuit", None],
}


def test_affinity_counts():
    result = affinity(COLLATION)
    ab = result["A", "B"]
    assert result["B", "A"] is ab
    assert (ab.comparable, ab.agreements, ab.shared, ab.differing) == (4, 3, 3.0, 1.0)
    assert ab.agreement == 0.75
    ad = result["A", "D"]
    assert (ad.comparable, ad.agreements, ad.shared, ad.differing) == (3, 1, 1.0, 2.0)
    # "in" is the reading of every witness (q = 1), so A and D have only
    # disagreements: the posterior is Beta(1, 9 + 2)
    assert ad.rho == pytest.approx(1 / 12, rel=1e-9)


def test_affinity_by_hand():
    # A and B: "in" (q = 1), "principio" (q = 1/3), "erat" (q = 2/3), one
    # disagreement; B and C: "sermo" has q = 1/2, as D has no reading there
    result = affinity(COLLATION)
    expected_ab = posterior([(1.0, 1.0), (1 / 3, 1.0), (2 / 3, 1.0)], 1.0)
    expected_bc = posterior([(1.0, 1.0), (2 / 3, 1.0), (1 / 2, 1.0)], 1.0)
    assert (
        result["A", "B"].rho,
        result["A", "B"].low,
        result["A", "B"].high,
    ) == expected_ab
    assert (
        result["B", "C"].rho,
        result["B", "C"].low,
        result["B", "C"].high,
    ) == expected_bc
    assert result["A", "B"].rho > result["B", "C"].rho


def test_affinity_counting_reading():
    # Disagreements count twice, agreements once
    result = affinity(COLLATION, counting="reading")
    assert (result["A", "B"].shared, result["A", "B"].differing) == (3.0, 2.0)
    assert (result["A", "D"].shared, result["A", "D"].differing) == (1.0, 4.0)


def test_affinity_weights():
    # Weight 0 for the reading of A (as a lemma): agreements in it are no
    # evidence, and disagreements with it count half in `locus` counting
    lemma = COLLATION["A"]
    result = affinity(COLLATION, weight=lambda i, r: 0.0 if r == lemma[i] else 1.0)
    ab = result["A", "B"]
    assert (ab.shared, ab.differing, ab.agreements) == (0.0, 0.5, 3)
    cd = result["C", "D"]
    # "erat" (weight 0) against "fuit" (weight 1)
    assert (cd.shared, cd.differing) == (1.0, 0.5)
    reading = affinity(
        COLLATION,
        weight=lambda i, r: 0.0 if r == lemma[i] else 1.0,
        counting="reading",
    )
    assert (reading["A", "B"].shared, reading["A", "B"].differing) == (0.0, 1.0)


def test_affinity_no_evidence():
    result = affinity({"A": ["x", None], "B": [None, "y"], "C": ["x", "y"]})
    assert result["A", "B"].comparable == 0
    assert math.isnan(result["A", "B"].rho)
    assert result["A", "B"].agreement is None
    # All weights zero
    result = affinity({"A": ["x"], "B": ["x"]}, weight=lambda i, r: 0.0)
    assert math.isnan(result["A", "B"].rho)
    assert result["A", "B"].agreement == 1.0


def test_affinity_witnesses_keep_frequencies():
    # Restricting the pairs does not change the frequencies of readings
    full = affinity(COLLATION)
    part = affinity(COLLATION, witnesses=["C", "B"])
    assert set(part) == {("B", "C"), ("C", "B")}
    assert part["B", "C"] == full["B", "C"]


def test_affinity_symmetric_and_omission_is_a_reading():
    coll = {"A": ["", "b", "c"], "B": ["", "b", "x"], "C": ["a", "b", "c"]}
    result = affinity(coll)
    assert result["A", "B"].agreements == 2
    assert result["A", "B"] == result["B", "A"]


@pytest.mark.parametrize(
    "coll,kwargs",
    [
        ({"A": ["x"]}, {}),
        ({"A": ["x"], "B": ["x", "y"]}, {}),
        (COLLATION, {"witnesses": ["A", "Z"]}),
        (COLLATION, {"counting": "pair"}),
        (COLLATION, {"prior": (0, 1)}),
        (COLLATION, {"weight": lambda i, r: -1.0}),
    ],
)
def test_affinity_errors(coll, kwargs):
    with pytest.raises(ValueError):
        affinity(coll, **kwargs)


def test_affinity_matrix():
    sigla, matrix = affinity_matrix(COLLATION, witnesses=["B", "A", "C"])
    assert sigla == ["B", "A", "C"]
    full = affinity(COLLATION)
    assert matrix[0][0] == 1.0
    assert matrix[0][1] == matrix[1][0] == full["A", "B"].rho
    assert matrix[2][1] == full["A", "C"].rho


def test_affinity_result_type():
    value = affinity(COLLATION)["A", "C"]
    assert isinstance(value, Affinity)
    assert 0.0 < value.low < value.rho < value.high < 1.0
    assert collation.MISSING is None


# Equivalence with CAB


def _cab_cases():
    data = json.loads(DATA.read_text(encoding="utf-8"))
    return data


def test_cab_posterior_reference():
    for case in _cab_cases()["posteriors"]:
        got = posterior([tuple(x) for x in case["shared"]], case["differing"], grid=301)
        assert got == pytest.approx(
            (case["mean"], case["low"], case["high"]), rel=1e-12, abs=1e-15
        )


@pytest.mark.parametrize("preset", ["var", "varreading", "all"])
def test_cab_pair_matrix_reference(preset):
    data = _cab_cases()
    types = data["types"]
    reading_weight = {
        "lem": 0.0,
        "var": 1.0,
        "reading": 0.0 if preset == "var" else 0.25,
    }
    result = affinity(
        data["collation"],
        weight=lambda i, r: reading_weight[types[i][r]],
        counting="reading",
        grid=301,
    )
    for pair, expected in data["pairs"][preset].items():
        a, b = pair.split("|")
        got = result[a, b]
        if expected["rho"] is None:
            assert math.isnan(got.rho)
            continue
        assert got.shared + got.differing == pytest.approx(expected["evidence"])
        assert (got.rho, got.low, got.high) == pytest.approx(
            (expected["rho"], expected["low"], expected["high"]), rel=1e-12, abs=1e-15
        )


def test_cab_transliteration_reference():
    # CAB's transliteration mode (pairMsStats): one observation per locus
    data = _cab_cases()["transliteration"]
    collation_ = {
        siglum: [None if code == 0 else code for code in codes]
        for siglum, codes in data["vectors"].items()
    }
    result = affinity(collation_, grid=301)
    for pair, expected in data["pairs"].items():
        a, b = pair.split("|")
        got = result[a, b]
        assert (got.comparable, got.agreements) == (
            expected["comparable"],
            expected["agreements"],
        )
        assert (got.rho, got.low, got.high) == pytest.approx(
            (expected["rho"], expected["low"], expected["high"]), rel=1e-12, abs=1e-15
        )


def test_cab_apparatus_subset_reference():
    # In apparatus mode, q is counted among all witnesses, also those not
    # selected
    data = _cab_cases()["apparatus_subset"]
    types = data["types"]
    reading_weight = {"lem": 0.0, "var": 1.0, "reading": 0.25}
    result = affinity(
        data["collation"],
        witnesses=data["selected"],
        weight=lambda i, r: reading_weight[types[i][r]],
        counting="reading",
        grid=301,
    )
    for pair, expected in data["pairs"].items():
        a, b = pair.split("|")
        got = result[a, b]
        assert (got.rho, got.low, got.high) == pytest.approx(
            (expected["rho"], expected["low"], expected["high"]), rel=1e-12, abs=1e-15
        )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"weight": lambda i, r: math.nan},
        {"weight": lambda i, r: math.inf},
        {"witnesses": ["A", "A"]},
    ],
)
def test_affinity_rejects_invalid(kwargs):
    with pytest.raises(ValueError):
        affinity(COLLATION, **kwargs)


@pytest.mark.parametrize("bad", [math.nan, math.inf])
def test_posterior_rejects_non_finite(bad):
    with pytest.raises(ValueError):
        posterior([], bad)
    with pytest.raises(ValueError):
        posterior([(0.5, bad)], 1.0)
