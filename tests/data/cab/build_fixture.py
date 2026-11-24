"""Convert the output of CAB's functions (make_fixture.sh) to the test fixture."""

import json
import sys
from pathlib import Path

WITNESSES = list("ABCDEFGH")
SHA = "b70c72091346770a904df9adec89817c0785dd1f91290099a3aaf7274ff98cc8"


def apparatus(loci, all_types=False):
    collation = {w: [None] * len(loci) for w in WITNESSES}
    types = []
    for i, locus in enumerate(loci):
        found = {}
        for value, kinds, witnesses in locus["r"]:
            if all_types:
                found.setdefault(value, set()).add(kinds[0])
            else:
                found[value] = kinds[0]
            for w in witnesses:
                collation[w][i] = value
        types.append({k: sorted(v) for k, v in found.items()} if all_types else found)
    return collation, types


def pairs(results, counts=False):
    out = {}
    for key, v in results.items():
        ci = v["bayesCI"]
        row = {
            "rho": v["bayes"],
            "low": ci[0] if ci else None,
            "high": ci[1] if ci else None,
        }
        if counts:
            row.update(comparable=v["comparable"], agreements=v["agreements"])
        else:
            row["evidence"] = v["bayesN"]
        out[key] = row
    return out


def main(work, target):
    work = Path(work)
    d1 = json.loads((work / "cab_out.json").read_text())
    d2 = json.loads((work / "cab_out2.json").read_text())
    d3 = json.loads((work / "cab_out3.json").read_text())
    collation, types = apparatus(d1["loci"])
    sub_collation, sub_types = apparatus(d2["loci"])
    all_collation, all_types = apparatus(d3["loci"], all_types=True)
    subset = {
        k: {kk: vv for kk, vv in v.items() if kk != "evidence"}
        for k, v in pairs(d2["app"]).items()
    }
    allp = {
        k: {kk: vv for kk, vv in v.items() if kk != "evidence"}
        for k, v in pairs(d3["app"]).items()
    }
    fixture = {
        "_comment": (
            "Synthetic apparatus (400 loci, 8 witnesses; readings r0 = lemma, others var or reading) and the values computed for it "
            "by the functions appPairMatrix and bayesPosterior of CAB (Corpus Avesticum Berolinense, Proximity according to shared "
            f"readings, common.js version 30 (loaded by the page of version 34; SHA-256 {SHA}), "
            "https://ada.geschkult.fu-berlin.de/manuscript_analysis/), run on 24 November 2026 with node and identity normalization. "
            "The sections 'transliteration' (pairMsStats, codes with 0 for missing) and 'apparatus_subset' (appPairMatrix with only the "
            "witnesses in 'selected', preset varreading) were computed the same way. The section 'apparatus_all' (preset all, with types "
            "var, om, reading, orth, and readings split into rows of different types that CAB merges into one group) was computed the same "
            "way; each reading lists the types of its rows. Rebuilt by tests/data/cab/make_fixture.sh."
        ),
        "collation": collation,
        "types": types,
        "pairs": {p: pairs(r) for p, r in d1["results"].items()},
        "posteriors": [
            {
                "shared": t["qs"],
                "differing": t["mm"],
                "mean": t["post"]["mean"],
                "low": t["post"]["lo"],
                "high": t["post"]["hi"],
            }
            for t in d1["tests"]
        ],
        "transliteration": {
            "vectors": d2["vectors"],
            "pairs": pairs(d2["ms"], counts=True),
        },
        "apparatus_subset": {
            "collation": sub_collation,
            "types": sub_types,
            "selected": d2["sub"],
            "pairs": subset,
        },
        "apparatus_all": {
            "collation": all_collation,
            "types": all_types,
            "pairs": allp,
        },
    }
    for section in ("transliteration",):
        for v in fixture[section]["pairs"].values():
            v.pop("evidence", None)
    Path(target).write_text(json.dumps(fixture, separators=(",", ":")))


if __name__ == "__main__":
    main(*sys.argv[1:])
