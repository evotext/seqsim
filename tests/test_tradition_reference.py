"""
test_tradition_reference
========================

Checks that `seqsim.tradition` and the order measures reproduce the reference
implementation of the Apophthegmata project (`reference_apophthegmata.py`)
exactly, on randomly generated traditions: a shared pool of items in a
conserved order, divided in chapters, with random block moves, losses,
repetitions, fragments with lacunae, and codices split into parts.

The tests using pandas (character matrices) or SciPy (Kendall's tau) are
skipped when these packages are not installed.
"""

# Import Python standard libraries
import math
import random

import pytest

# Import the library being tested
from seqsim import order, token, tradition

import reference_apophthegmata as ref

SEEDS = range(12)


def make_tradition(seed):
    """
    Returns witnesses, layouts, and Monastica-like codes of a random tradition.
    """

    rng = random.Random(seed)
    n_items = rng.randint(40, 120)
    n_chapters = rng.randint(3, 8)
    pool = [f"s{idx:03d}" for idx in range(n_items)]
    chapter = {
        item: f"C{idx * n_chapters // n_items + 1:02d}" for idx, item in enumerate(pool)
    }

    def evolve(seq):
        seq = list(seq)
        for _ in range(rng.randint(0, 3)):  # block moves
            if len(seq) < 4:
                break
            i, j = sorted(rng.sample(range(len(seq)), 2))
            block, rest = seq[i:j], seq[:i] + seq[j:]
            k = rng.randint(0, len(rest))
            seq = rest[:k] + block + rest[k:]
        seq = [item for item in seq if rng.random() > rng.choice([0.0, 0.05, 0.3])]
        if seq and rng.random() < 0.3:  # a repeated item
            seq.insert(rng.randint(0, len(seq)), rng.choice(seq))
        if rng.random() < 0.4:  # a fragment
            i = rng.randint(0, len(seq) // 3)
            j = rng.randint(2 * len(seq) // 3, len(seq))
            seq = seq[i:j]
        if rng.random() < 0.2:  # an excerpt collection
            seq = sorted(rng.sample(seq, min(len(seq), 8)), key=pool.index)
        return seq

    witnesses, layouts, codes = [], {}, {}
    n_codices = rng.randint(6, 14)
    for codex_idx in range(n_codices):
        codex = f"Codex {codex_idx}"
        stories = evolve(pool)
        if not stories:
            continue
        kind = rng.random()
        if kind < 0.6 or len(stories) < 6:
            parts = [(None, stories)]
        elif kind < 0.8:  # consecutive parts
            cut = rng.randint(1, len(stories) - 1)
            parts = [("A", stories[:cut]), ("B", stories[cut:])]
        else:  # interleaved parts, with a second copy of a section
            parts = [("A", stories[0::2]), ("C", stories[1::2])]
            parts.append(("B", stories[: max(2, len(stories) // 3)]))
        code = rng.choice(["GS", "GS", "PJ"])
        for part, part_stories in parts:
            wid = codex if part is None else f"{codex} {part}"
            witnesses.append(ref.Witness(wid, list(part_stories), codex, part))
            codes[wid] = code if rng.random() > 0.1 else "PJ"
            layouts[wid] = make_layout(rng, part_stories, chapter)

    return witnesses, layouts, codes


def make_layout(rng, stories, chapter):
    """
    Returns a Monastica-like layout, in paragraphs, of a witness.
    """

    layout = {}
    if rng.random() < 0.15:
        return layout
    n = len(stories)
    paragraphs = max(1, round(n / rng.choice([1.0, 1.2, 1.5])))
    layout["paragraphs"] = paragraphs
    events = []
    if rng.random() < 0.85:
        previous = None
        for idx, item in enumerate(stories):
            if chapter[item] != previous and (
                previous is None or chapter[item] > previous
            ):
                pos = max(0, round(idx * paragraphs / n) + rng.choice([-1, 0, 0, 1]))
                events.append(
                    {"kind": "chapter", "paragraph": pos, "label": chapter[item]}
                )
                previous = chapter[item]
    for _ in range(rng.choice([0, 0, 1, 2])):
        events.append({"kind": "lacuna", "paragraph": rng.randint(0, paragraphs)})
    rng.shuffle(events)  # the order of events of different kinds is irrelevant
    events.sort(key=lambda event: event["kind"] != "chapter")
    layout["events"] = events
    return layout


def voters_of(witnesses, layouts):
    return [
        (
            w.stories,
            int(layouts.get(w.id, {}).get("paragraphs", 0)),
            [
                (e["paragraph"], e["label"])
                for e in ref._events(layouts.get(w.id, {}), "chapter")
            ],
        )
        for w in witnesses
    ]


def gaps_of(witness, layout):
    paragraphs = int(layout.get("paragraphs", 0))
    lacunae = [e["paragraph"] for e in ref._events(layout, "lacuna")]
    return tradition.proportional_positions(len(witness.order), paragraphs, lacunae)


def build(seed):
    """
    Runs both implementations of the frame and coverage on a tradition.
    """

    witnesses, layouts, codes = make_tradition(seed)
    voters = [w for w in witnesses if codes[w.id] == "GS"]
    frame_witnesses = [w for w in witnesses if not w.id.endswith(" B")]

    reference, support = ref.reference_chapters(voters, layouts)
    label_of, new_support = tradition.reference_labels(voters_of(voters, layouts))
    assert list(label_of.items()) == list(reference.items())
    assert list(new_support.items()) == list(support.items())

    frame = ref.consensus_order(frame_witnesses)
    assert tradition.consensus_order([w.stories for w in frame_witnesses]) == frame

    chapter_of = ref.fill_reference(frame, reference)
    filled = tradition.fill_labels(frame, label_of)
    assert list(filled.items()) == list(chapter_of.items())

    return witnesses, layouts, codes, frame, chapter_of


@pytest.mark.parametrize("seed", SEEDS)
def test_labels(seed):
    witnesses, layouts, _, frame, _ = build(seed)
    for w in witnesses:
        layout = layouts[w.id]
        paragraphs = int(layout.get("paragraphs", 0))
        events = [(e["paragraph"], e["label"]) for e in ref._events(layout, "chapter")]
        prior = ref.proportional_labels(len(w.stories), layout)
        assert (
            tradition.proportional_labels(len(w.stories), paragraphs, events) == prior
        )
        assert gaps_of(w, layout) == ref.lacuna_story_indices(len(w.order), layout)

        # Monotone labels against a noisy reference
        rng = random.Random(w.id)
        noisy = {
            item: rng.choice(["C01", "C02", "C03"])
            for item in frame
            if rng.random() < 0.5
        }
        label_order = [label for _, label in events]
        assert tradition.monotone_labels(
            w.stories, label_order, noisy, prior
        ) == ref.fit_labels(w.stories, label_order, noisy, prior)
        assert tradition.monotone_labels(
            w.stories, label_order, noisy, prior, prior_weight=1.5
        ) == ref.fit_labels(w.stories, label_order, noisy, prior, prior_weight=1.5)


@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize("min_share", [0.06, 0.3])
def test_coverage(seed, min_share):
    witnesses, layouts, codes, frame, chapter_of = build(seed)
    for w in witnesses:
        layout = layouts[w.id]
        expected = ref.coverage(
            w, frame, chapter_of, layout or None, min_share=min_share
        )
        cells = tradition.coverage(
            w.stories, frame, chapter_of, gaps=gaps_of(w, layout), min_share=min_share
        )
        assert list(cells.items()) == list(expected.items())
        assert tradition.covered_labels(
            w.stories, frame, chapter_of, min_share=min_share
        ) == ref.covered_chapters(w, frame, chapter_of, min_share=min_share)
        assert tradition.density(cells) == ref.density(list(expected.values()))


def merged_manuscripts(witnesses, codes, max_overlap=0.2):
    """
    Merges the parts of each codex with `merge_parts`, as the project does.
    """

    by_codex = {}
    for w in witnesses:
        by_codex.setdefault(w.codex, []).append(w)
    result, parts_of = [], {}
    for codex, parts in by_codex.items():
        by_label = {w.part: w for w in parts}
        merged, separate = tradition.merge_parts(
            [(w.part, w.stories) for w in parts],
            max_overlap=max_overlap,
            compatible=lambda a, b: codes.get(by_label[a].id)
            == codes.get(by_label[b].id),
        )
        if merged:
            label = "+".join(part for part, _ in merged)
            stories = [item for _, items in merged for item in items]
            witness = ref.Witness(f"{codex} {label}", stories, codex, label)
            parts_of[witness.id] = [by_label[part] for part, _ in merged]
            result.append(witness)
        result.extend(by_label[part] for part, _ in separate)
    return result, parts_of


@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize("max_overlap", [0.2, 0.6])
def test_manuscripts(seed, max_overlap):
    witnesses, layouts, codes, frame, chapter_of = build(seed)
    expected, expected_parts = ref.core_manuscripts(witnesses, codes, max_overlap)
    result, parts_of = merged_manuscripts(witnesses, codes, max_overlap)

    assert [(w.id, w.stories) for w in result] == [(w.id, w.stories) for w in expected]
    assert {k: [w.id for w in v] for k, v in parts_of.items()} == {
        k: [w.id for w in v] for k, v in expected_parts.items()
    }

    for w in result:
        parts = parts_of.get(w.id, [w])
        cells = tradition.merge_coverage(
            [
                tradition.coverage(
                    p.stories, frame, chapter_of, gaps=gaps_of(p, layouts.get(p.id, {}))
                )
                for p in parts
            ],
            frame,
        )
        expected_cells = ref.core_coverage(
            w, expected_parts, frame, chapter_of, layouts
        )
        assert list(cells.items()) == list(expected_cells.items())


def matrices(seed):
    witnesses, layouts, codes, frame, chapter_of = build(seed)
    taxa, parts_of = merged_manuscripts(witnesses, codes)
    taxa = [w for w in taxa if len(w) >= 5]
    coverages = {}
    for w in taxa:
        parts = parts_of.get(w.id, [w])
        coverages[w.id] = tradition.merge_coverage(
            [
                tradition.coverage(
                    p.stories, frame, chapter_of, gaps=gaps_of(p, layouts.get(p.id, {}))
                )
                for p in parts
            ],
            frame,
        )
    return taxa, parts_of, frame, chapter_of, coverages


@pytest.mark.parametrize("seed", SEEDS)
def test_characters(seed):
    pd = pytest.importorskip("pandas")
    taxa, parts_of, frame, chapter_of, coverages = matrices(seed)

    content = tradition.content_characters(coverages)
    expected_content = ref.content_characters(coverages)
    assert list(content.characters) == list(expected_content.columns)
    assert list(content.taxa) == list(expected_content.index)
    assert tradition.to_phylip(content) == ref.to_phylip(expected_content)

    witnesses = {w.id: w.stories for w in taxa}
    parts = {wid: [p.stories for p in ps] for wid, ps in parts_of.items()}
    adjacency = tradition.adjacency_characters(witnesses, frame, parts_of=parts)
    expected_adjacency = ref.adjacency_characters(taxa, frame, parts_of=parts_of)
    assert list(adjacency.characters) == list(expected_adjacency.columns)
    assert list(adjacency.taxa) == list(expected_adjacency.index)
    for name in adjacency.characters:
        assert adjacency.column(name) == list(expected_adjacency[name])
    assert tradition.to_phylip(adjacency) == ref.to_phylip(expected_adjacency)

    # The conversion to pandas gives the same data frame
    converted = pd.DataFrame(
        adjacency.as_dict(), index=list(adjacency.taxa), dtype=object
    )
    pd.testing.assert_frame_equal(converted, expected_adjacency)

    # Combined matrix and partitions
    combined = tradition.concat_characters([content, adjacency])
    expected_combined = pd.concat([expected_content, expected_adjacency], axis=1)
    assert tradition.to_phylip(combined) == ref.to_phylip(expected_combined)
    charsets = tradition.charset_ranges(
        {"content": content.shape[1], "adjacency": adjacency.shape[1]}
    )
    assert tradition.partition_nexus(charsets) == ref.partitions_text(
        expected_content, expected_adjacency
    )
    assert tradition.to_nexus(combined, charsets) == ref.nexus(
        expected_combined, charsets
    )
    assert tradition.to_nexus(combined) == ref.nexus(expected_combined)

    # Restriction to some taxa
    rng = random.Random(seed)
    rows = rng.sample(list(content.taxa), max(1, len(content.taxa) * 2 // 3))
    for matrix, expected in [
        (content, expected_content),
        (adjacency, expected_adjacency),
    ]:
        restricted = tradition.restrict(matrix, rows)
        expected_restricted = ref.restricted(expected, rows)
        assert list(restricted.characters) == list(expected_restricted.columns)
        assert tradition.to_phylip(restricted) == ref.to_phylip(expected_restricted)

    # Blocks for resampling
    position = {item: pos for pos, item in enumerate(frame)}
    for matrix, kind in [(content, "content"), (adjacency, "adjacency")]:
        expected_blocks = ref.character_blocks(
            list(matrix.characters), position, chapter_of, kind
        )
        separator = ">" if kind == "adjacency" else None
        assert tradition.character_blocks(
            matrix.characters, position, size=20, separator=separator
        ) == list(expected_blocks["block20"])
        assert tradition.character_blocks(
            matrix.characters, label_of=chapter_of, separator=separator
        ) == list(expected_blocks["chapter"])


@pytest.mark.parametrize("seed", SEEDS)
def test_pairwise_measures(seed):
    witnesses, *_ = build(seed)
    for a in witnesses:
        for b in witnesses:
            assert token.containment(a.order, b.order) == ref.containment(a, b)
            ra, rb = order.restrict_to_shared(a.order, b.order)
            assert (ra, rb) == ref._shared_orders(a, b)
            assert (ra, rb) == order.restrict_to_shared(
                a.stories, b.stories, repeats="first"
            )
            agreement = ref.adjacency_agreement(a, b, min_shared=2)
            if math.isnan(agreement):
                assert len(ra) < 2
            else:
                assert order.breakpoint_simil(ra, rb, boundaries=False) == agreement


@pytest.mark.parametrize("seed", SEEDS)
def test_kendall_tau_against_scipy(seed):
    pytest.importorskip("scipy")
    witnesses, *_ = build(seed)
    for a in witnesses:
        for b in witnesses:
            expected = ref.order_tau(a, b, min_shared=2)
            value = order.kendall_tau_simil(a.order, b.order)
            if math.isnan(expected):
                assert math.isnan(value)
            else:
                assert value == pytest.approx(expected, abs=1e-12)


def test_mrbayes_lines():
    pd = pytest.importorskip("pandas")
    dates = pd.DataFrame(
        {
            "taxon": ["W_1", "W_2", "W_3"],
            "age_min": [500, 300, 120],
            "age_max": [500, 400, 150],
        }
    )
    tip_ages = {r.taxon: (r.age_min, r.age_max) for _, r in dates.iterrows()}
    assert tradition.mrbayes_calibrations(tip_ages) == ref.calibrate_lines(dates)

    groups = {
        "west": (["W_1", "W_2"], (600, 900)),
        "east": (["W_3", "W_4"], (700, 800)),
    }
    assert tradition.mrbayes_constraints(groups) == ref.constraint_lines(groups)
    assert tradition.mrbayes_constraints({}) == ref.constraint_lines({}) == []
