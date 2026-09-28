"""
Reference implementation from the Apophthegmata project, for equivalence tests.

The functions below are copied from `tresoldi/apophthegmata` at commit
`65dd2c8` (`apophthegmata/corpus.py`, `apophthegmata/coverage.py`,
`apophthegmata/characters.py`, and helpers of the analysis scripts), with
only the changes needed to run them outside the project:

  * `Witness` is a minimal stand-in for the project's class: `content` is the
    set of its stories, `order` their first occurrences, and `len()` the
    number of distinct stories, as in the project (`apophthegmata/corpus.py`);
  * `SystematicCore.coverage` and `SystematicCore.manuscripts` are reduced
    to the functions `core_coverage` and `core_manuscripts`, taking the
    frame, the labels, the layouts, and the codes as arguments;
  * the pandas and SciPy functions need those packages, imported only when
    available; the tests using them are skipped otherwise.

Do not "fix" this file: it documents the behaviour that `seqsim.tradition`
and the order measures must reproduce.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    pd = None

try:
    import numpy as np
    from scipy.stats import kendalltau
except ImportError:  # pragma: no cover
    np = None
    kendalltau = None


@dataclass
class Witness:
    id: str
    stories: list = field(default_factory=list)
    codex: str = ""
    part: str | None = None

    @property
    def content(self) -> set:
        return set(self.stories)

    @property
    def order(self) -> list:
        return list(dict.fromkeys(self.stories))

    def __len__(self) -> int:
        return len(self.content)


# ---------------------------------------------------------------------------
# apophthegmata/corpus.py


def containment(a: Witness, b: Witness) -> float:
    return len(a.content & b.content) / len(a.content) if a.content else 0.0


def _shared_orders(a: Witness, b: Witness) -> tuple[list[str], list[str]]:
    shared = a.content & b.content
    return [s for s in a.order if s in shared], [s for s in b.order if s in shared]


def order_tau(a: Witness, b: Witness, min_shared: int = 20) -> float:
    ra, rb = _shared_orders(a, b)
    if len(ra) < min_shared:
        return float("nan")
    rank_b = {s: i for i, s in enumerate(rb)}
    return float(kendalltau(np.arange(len(ra)), [rank_b[s] for s in ra])[0])


def adjacency_agreement(a: Witness, b: Witness, min_shared: int = 20) -> float:
    ra, rb = _shared_orders(a, b)
    if len(ra) < min_shared:
        return float("nan")
    adjacencies_a = set(zip(ra, ra[1:], strict=False))
    adjacencies_b = set(zip(rb, rb[1:], strict=False))
    return len(adjacencies_a & adjacencies_b) / (len(ra) - 1)


# ---------------------------------------------------------------------------
# apophthegmata/coverage.py

Layout = Mapping[str, Any]


def _events(layout: Layout, kind: str) -> list[dict]:
    return [e for e in layout.get("events", []) if e["kind"] == kind]


def proportional_labels(n_stories: int, layout: Layout) -> list[str | None]:
    paragraphs = int(layout.get("paragraphs", 0))
    chapters = _events(layout, "chapter")
    if not paragraphs or not chapters:
        return [None] * n_stories
    labels: list[str | None] = []
    for i in range(n_stories):
        p = i * paragraphs / n_stories
        label = None
        for event in chapters:
            if event["paragraph"] <= p:
                label = event["label"]
        labels.append(label)
    return labels


def fit_labels(
    stories: Sequence[str],
    chapter_order: Sequence[str],
    reference: Mapping[str, str],
    prior: Sequence[str | None],
    prior_weight: float = 0.25,
) -> list[str]:
    n, k = len(stories), len(chapter_order)
    if not n or not k:
        return [chapter_order[0] if k else ""] * n

    def gain(i: int, j: int) -> float:
        chapter = chapter_order[j]
        score = 1.0 if reference.get(stories[i]) == chapter else 0.0
        return score + (prior_weight if prior[i] == chapter else 0.0)

    best = [[0.0] * k for _ in range(n)]
    back = [[0] * k for _ in range(n)]
    for j in range(k):
        best[0][j] = gain(0, j)
    for i in range(1, n):
        running, arg = best[i - 1][0], 0
        for j in range(k):
            if best[i - 1][j] > running:
                running, arg = best[i - 1][j], j
            best[i][j] = running + gain(i, j)
            back[i][j] = arg
    j = max(range(k), key=lambda c: best[n - 1][c])
    labels = [""] * n
    for i in range(n - 1, -1, -1):
        labels[i] = chapter_order[j]
        j = back[i][j]
    return labels


def reference_chapters(
    witnesses: Iterable[Witness],
    layouts: Mapping[str, Layout],
    rounds: int = 3,
) -> tuple[dict[str, str], dict[str, float]]:
    voters = [w for w in witnesses if _events(layouts.get(w.id, {}), "chapter")]
    labels = {w.id: proportional_labels(len(w.stories), layouts[w.id]) for w in voters}
    reference: dict[str, str] = {}
    support: dict[str, float] = {}
    for _ in range(rounds):
        votes: dict[str, Counter[str]] = defaultdict(Counter)
        for w in voters:
            for story, label in zip(w.stories, labels[w.id], strict=True):
                if label:
                    votes[story][label] += 1
        reference = {s: v.most_common(1)[0][0] for s, v in votes.items()}
        support = {
            s: v.most_common(1)[0][1] / sum(v.values()) for s, v in votes.items()
        }
        for w in voters:
            order = [e["label"] for e in _events(layouts[w.id], "chapter")]
            prior = proportional_labels(len(w.stories), layouts[w.id])
            labels[w.id] = list(fit_labels(w.stories, order, reference, prior))
    return reference, support


def consensus_order(witnesses: Sequence[Witness]) -> list[str]:
    frame: list[str] = []
    placed: set[str] = set()
    for witness in witnesses:
        anchor = -1
        insertions: dict[int, list[str]] = defaultdict(list)
        index = {s: i for i, s in enumerate(frame)}
        for story in witness.order:
            if story in placed:
                anchor = index[story]
            else:
                insertions[anchor].append(story)
                placed.add(story)
        merged = list(insertions.get(-1, []))
        for i, story in enumerate(frame):
            merged.append(story)
            merged.extend(insertions.get(i, []))
        frame = merged
    return frame


def fill_reference(
    frame: Sequence[str], reference: Mapping[str, str]
) -> dict[str, str]:
    filled = dict(reference)
    last = None
    pending: list[str] = []
    for story in frame:
        if story in reference:
            last = reference[story]
            for p in pending:
                filled[p] = last
            pending = []
        elif last is not None:
            filled[story] = last
        else:
            pending.append(story)
    return filled


def lacuna_story_indices(n_stories: int, layout: Layout) -> list[int]:
    paragraphs = int(layout.get("paragraphs", 0))
    if not paragraphs:
        return []
    return [
        round(e["paragraph"] * n_stories / paragraphs)
        for e in _events(layout, "lacuna")
    ]


def covered_chapters(
    witness: Witness,
    frame: Sequence[str],
    chapter_of: Mapping[str, str],
    min_stories: int = 2,
    min_share: float = 0.06,
) -> set[str]:
    chapter_size = Counter(chapter_of.get(s) for s in frame)
    in_frame = set(frame)
    per_chapter = Counter(chapter_of.get(s) for s in witness.content if s in in_frame)
    return {
        c
        for c, count in per_chapter.items()
        if c is not None
        and count >= min_stories
        and count >= min_share * chapter_size[c]
    }


def coverage(
    witness: Witness,
    frame: Sequence[str],
    chapter_of: Mapping[str, str],
    layout: Layout | None = None,
    min_stories: int = 2,
    min_share: float = 0.06,
) -> dict[str, int | None]:
    position = {s: i for i, s in enumerate(frame)}
    present = [s for s in witness.order if s in position]
    if not present:
        return dict.fromkeys(frame)

    covered = covered_chapters(witness, frame, chapter_of, min_stories, min_share)
    first = min(position[s] for s in present)
    last = max(position[s] for s in present)

    lost: set[int] = set()
    if layout:
        sequence = witness.order
        for cut in lacuna_story_indices(len(sequence), layout):
            before = [position[s] for s in sequence[:cut] if s in position]
            after = [position[s] for s in sequence[cut:] if s in position]
            lo = before[-1] if before else -1
            hi = after[0] if after else len(frame)
            if lo < hi:
                lost.update(range(lo + 1, hi))

    content = witness.content
    result: dict[str, int | None] = {}
    for i, story in enumerate(frame):
        if story in content:
            result[story] = 1
        elif chapter_of.get(story) not in covered or i < first or i > last or i in lost:
            result[story] = None
        else:
            result[story] = 0
    return result


def core_coverage(
    witness: Witness,
    parts_of: Mapping[str, list[Witness]],
    frame: Sequence[str],
    chapter_of: Mapping[str, str],
    layouts: Mapping[str, Layout],
    **kwargs: Any,
) -> dict[str, int | None]:
    """`SystematicCore.coverage`."""
    parts = parts_of.get(witness.id, [witness])
    per_part = [
        coverage(part, frame, chapter_of, layouts.get(part.id) or None, **kwargs)
        for part in parts
    ]
    if len(per_part) == 1:
        return per_part[0]
    combined: dict[str, int | None] = {}
    for story in frame:
        values = [cells[story] for cells in per_part]
        combined[story] = 1 if 1 in values else 0 if 0 in values else None
    return combined


def core_manuscripts(
    witnesses: Sequence[Witness], codes: Mapping[str, str], max_overlap: float = 0.2
) -> tuple[list[Witness], dict[str, list[Witness]]]:
    """`SystematicCore.manuscripts`, without the final sort and metadata."""
    parts_of: dict[str, list[Witness]] = {}
    by_codex: dict[str, list[Witness]] = defaultdict(list)
    for w in witnesses:
        by_codex[w.codex].append(w)
    result = []
    for codex, parts in by_codex.items():
        if len(parts) == 1:
            result.append(parts[0])
            continue
        merged = [max(parts, key=len)]
        content = set(merged[0].content)
        separate = []
        for part in sorted(parts, key=len, reverse=True)[1:]:
            overlap = len(part.content & content) / len(part.content)
            if overlap < max_overlap and codes.get(part.id) == codes.get(merged[0].id):
                merged.append(part)
                content |= part.content
            else:
                separate.append(part)
        if len(merged) == 1:
            result.extend(parts)
            continue
        merged.sort(key=lambda w: w.part or "")
        label = "+".join(w.part or "" for w in merged)
        witness = Witness(
            id=f"{codex} {label}",
            stories=[s for w in merged for s in w.stories],
            codex=codex,
            part=label,
        )
        parts_of[witness.id] = merged
        result.append(witness)
        result.extend(separate)
    return result, parts_of


# ---------------------------------------------------------------------------
# apophthegmata/characters.py (pandas)

Cell = int | None


def informative(column: Iterable[Cell], min_each: int = 2) -> bool:
    values = list(column)
    return values.count(1) >= min_each and values.count(0) >= min_each


def content_characters(coverage: Mapping[str, Mapping[str, Cell]], min_each: int = 2):
    frame = pd.DataFrame(coverage).T
    keep = [c for c in frame.columns if informative(frame[c], min_each)]
    return frame[keep]


def adjacencies(witness: Witness, universe: set[str]) -> set[tuple[str, str]]:
    order = [s for s in witness.order if s in universe]
    return set(zip(order, order[1:], strict=False))


def adjacency_characters(
    witnesses: Sequence[Witness],
    universe: Sequence[str],
    min_each: int = 2,
    parts_of: Mapping[str, Sequence[Witness]] | None = None,
):
    parts_of = parts_of or {}
    universe_set = set(universe)
    parts = {w.id: list(parts_of.get(w.id, [w])) for w in witnesses}
    per_witness = {
        wid: set().union(*(adjacencies(p, universe_set) for p in ps))
        for wid, ps in parts.items()
    }
    candidates: dict[tuple[str, str], int] = {}
    for adj in per_witness.values():
        for pair in adj:
            candidates[pair] = candidates.get(pair, 0) + 1
    columns = {}
    for (x, y), count in sorted(candidates.items()):
        if count < min_each:
            continue
        cells: list[Cell] = []
        for w in witnesses:
            if (x, y) in per_witness[w.id]:
                cells.append(1)
            elif any(x in p.content and y in p.content for p in parts[w.id]):
                cells.append(0)
            else:
                cells.append(None)
        if informative(cells, min_each):
            columns[f"{x}>{y}"] = cells
    return pd.DataFrame(columns, index=[w.id for w in witnesses], dtype=object)


def taxon_name(witness_id: str) -> str:
    return witness_id.replace(" ", "_").replace("+", "-")


def to_phylip(matrix) -> str:
    symbol = {1: "1", 0: "0"}
    lines = [f"{matrix.shape[0]} {matrix.shape[1]}"]
    for witness, row in matrix.iterrows():
        lines.append(
            f"{taxon_name(str(witness))} " + "".join(symbol.get(v, "?") for v in row)
        )
    return "\n".join(lines) + "\n"


def partitions_text(content, adjacency) -> str:
    """The partition file of `write_matrices` and `write_partitioned`."""
    combined = pd.concat([content, adjacency], axis=1)
    return (
        f"#nexus\nbegin sets;\n  charset content = 1-{content.shape[1]};\n"
        f"  charset adjacency = {content.shape[1] + 1}-{combined.shape[1]};\nend;\n"
    )


# ---------------------------------------------------------------------------
# Helpers of the analysis scripts

BLOCK = 20


def character_blocks(
    columns: list[str], position: dict[str, int], chapter: dict[str, str], kind: str
):
    """analysis/05_robustness.py"""
    first = [c.split(">")[0] if kind == "adjacency" else c for c in columns]
    return {
        "block20": np.array([position[s] // BLOCK for s in first]),
        "chapter": np.array([chapter[s] for s in first]),
    }


def density(values: list) -> float:
    """analysis/06_backbone_placement.py, `densities`, for one witness."""
    present, absent = values.count(1), values.count(0)
    return present / (present + absent) if present + absent else 0.0


def restricted(matrix, rows: list[str]):
    """analysis/07_prepare_bayesian.py (and `keep_informative` of 06)."""
    sub = matrix.loc[rows]
    return sub[[c for c in sub.columns if informative(sub[c])]]


def nexus(matrix, charsets: dict[str, tuple[int, int]] | None = None) -> str:
    """analysis/07_prepare_bayesian.py"""
    symbol = {1: "1", 0: "0"}
    names = [taxon_name(str(w)) for w in matrix.index]
    width = max(len(n) for n in names) + 2
    lines = [
        "#NEXUS",
        "begin data;",
        f"  dimensions ntax={matrix.shape[0]} nchar={matrix.shape[1]};",
        "  format datatype=restriction missing=? gap=-;",
        "  matrix",
    ]
    for name, (_, row) in zip(names, matrix.iterrows(), strict=True):
        lines.append(
            f"    {name.ljust(width)}" + "".join(symbol.get(v, "?") for v in row)
        )
    lines += ["  ;", "end;"]
    if charsets:
        lines += ["begin sets;"]
        lines += [f"  charset {name} = {a}-{b};" for name, (a, b) in charsets.items()]
        lines += ["end;"]
    return "\n".join(lines) + "\n"


def calibrate_lines(dates) -> list[str]:
    """analysis/07_prepare_bayesian.py"""
    lines = []
    for _, r in dates.iterrows():
        if r.age_min == r.age_max:
            lines.append(f"  calibrate {r.taxon} = fixed({r.age_min});")
        else:
            lines.append(f"  calibrate {r.taxon} = uniform({r.age_min},{r.age_max});")
    return lines


def constraint_lines(groups) -> list[str]:
    """The group lines of `dated_block` in analysis/07_prepare_bayesian.py."""
    lines = []
    if groups:
        for name, (members, _) in groups.items():
            lines.append(f"  constraint {name} = {' '.join(members)};")
        for name, (_, (young, old)) in groups.items():
            lines.append(f"  calibrate {name} = uniform({young},{old});")
        lines.append(f"  prset topologypr=constraints({','.join(groups)});")
    return lines
