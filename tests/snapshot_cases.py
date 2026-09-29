"""
Cases for the snapshot test of all measures (see `test_snapshot.py`).

The snapshot records the exact output of every public measure on these
inputs, so that refactoring cannot silently change any result. Regenerate
it only for intended changes of results, with:

    python tests/snapshot_cases.py
"""

# Import Python standard libraries
import json
import pathlib

# Import the library being tested
from seqsim import alignment, compression, edit, order, sequence, token

SNAPSHOT = pathlib.Path(__file__).parent / "data" / "snapshot.json"

# Pairs of sequences: strings, lists and tuples, mixed element types,
# repeated elements, disjoint and empty sequences, and a longer pair
PAIRS = [
    ("kitten", "sitting"),
    ("abcdef", "badcfe"),
    ("ca", "abc"),
    ("abacada", "acabada"),
    ("aab", "abb"),
    ("abcdeXXXXXfghij", "abcdefghij"),
    ("abcdefgh", "efghabcd"),
    ((1, 2, 3, 4, 5), [1, 2, 4, 3, 6, 7]),
    ((1, 2, 3), ["a", "b", "c", "d"]),
    ([1, (2, 3), None, "x"], [None, 1, "x", (2, 3), (2, 3)]),
    ("abc", "xyz"),
    ("a", ""),
    ("", ""),
    ("abc", "abc"),
    (
        list("the quick brown fox jumps over the lazy dog"),
        list("the quick brown dog jumps over the lazy fox"),
    ),
    ([i % 7 for i in range(40)], [(i * 3) % 7 for i in range(35)]),
]

# Functions and keyword arguments to record; `normal` is added for all
# functions that accept it
FUNCTIONS = {
    "alignment.nw_dissim": (alignment.nw_dissim, [{}, {"gap_open": 2.0}]),
    "alignment.sw_simil": (alignment.sw_simil, [{}, {"gap_open": 1.0}]),
    "compression.entropy_ncd_dissim": (compression.entropy_ncd_dissim, [{}]),
    "compression.lz76_dissim": (compression.lz76_dissim, [{}]),
    "compression.lzma_ncd_dissim": (compression.lzma_ncd_dissim, [{}]),
    "edit.birnbaum_dissim": (edit.birnbaum_dissim, [{}]),
    "edit.birnbaum_simil": (edit.birnbaum_simil, [{}]),
    "edit.block_move_dissim": (
        edit.block_move_dissim,
        [{}, {"directional": True}],
    ),
    "edit.bulk_delete_dist": (edit.bulk_delete_dist, [{}, {"max_del_len": 2}]),
    "edit.damerau_dist": (edit.damerau_dist, [{}]),
    "edit.damerau_gld_dist": (edit.damerau_gld_dist, [{}]),
    "edit.fragile_ends_dissim": (
        edit.fragile_ends_dissim,
        [{}, {"frag_start": 25.0, "frag_end": 0.0}],
    ),
    "edit.gst_dissim": (edit.gst_dissim, [{}, {"min_match": 1}]),
    "edit.indel_dist": (edit.indel_dist, [{}]),
    "edit.indel_gld_dist": (edit.indel_gld_dist, [{}]),
    "edit.jaro_dissim": (edit.jaro_dissim, [{}]),
    "edit.jaro_winkler_dissim": (edit.jaro_winkler_dissim, [{}]),
    "edit.lcs_dist": (edit.lcs_dist, [{}]),
    "edit.levenshtein_dist": (edit.levenshtein_dist, [{}]),
    "edit.levenshtein_gld_dist": (edit.levenshtein_gld_dist, [{}]),
    "edit.levenshtein_ned_dist": (edit.levenshtein_ned_dist, [{}]),
    "edit.mmcwpa_dissim": (edit.mmcwpa_dissim, [{}]),
    "edit.osa_dissim": (edit.osa_dissim, [{}]),
    "edit.stemmatological_dissim": (
        edit.stemmatological_dissim,
        [{}, {"frag_start": 20.0, "max_del_len": 2}],
    ),
    "order.block_interchange_dissim": (order.block_interchange_dissim, [{}]),
    "order.breakpoint_dissim": (
        order.breakpoint_dissim,
        [{}, {"boundaries": False}],
    ),
    "order.breakpoint_simil": (
        order.breakpoint_simil,
        [{}, {"boundaries": False}],
    ),
    "order.cayley_dissim": (order.cayley_dissim, [{}]),
    "order.footrule_dissim": (order.footrule_dissim, [{}, {"ell": 50}]),
    "order.iebp_estimate": (order.iebp_estimate, [{}, {"boundaries": False}]),
    "order.kendall_tau_dissim": (order.kendall_tau_dissim, [{}, {"p": 1.0}]),
    "order.kendall_tau_simil": (order.kendall_tau_simil, [{}]),
    "order.ulam_dist": (order.ulam_dist, [{}]),
    "sequence.ratcliff_obershelp_dissim": (
        sequence.ratcliff_obershelp_dissim,
        [{}],
    ),
    "token.containment": (token.containment, [{}, {"size": 2}]),
    "token.jaccard_dissim": (token.jaccard_dissim, [{}]),
    "token.qgram_dissim": (token.qgram_dissim, [{}, {"q": 3, "pad": False}]),
    "token.sorensen_dissim": (token.sorensen_dissim, [{}]),
    "token.subseq_jaccard_dissim": (token.subseq_jaccard_dissim, [{}]),
    "token.tversky_simil": (
        token.tversky_simil,
        [{}, {"alpha": 1.0, "beta": 0.0}],
    ),
}

# Sequences of sequences for Monge-Elkan
TITLE_PAIRS = [
    (["vita antonii", "passio agnetis"], ["passio s. agnetis", "vita antonij"]),
    (["a"], []),
    ([], []),
]


def _call(func, seq_x, seq_y, kwargs):
    try:
        return repr(func(seq_x, seq_y, **kwargs))
    except Exception as exception:  # recorded, so that errors are preserved too
        return f"{type(exception).__name__}"


def compute():
    """
    Returns the outputs of all cases, as a mapping from case keys to reprs.
    """

    import inspect

    results = {}
    for name, (func, variants) in FUNCTIONS.items():
        accepts_normal = "normal" in inspect.signature(func).parameters
        for kwargs in variants:
            options = [kwargs]
            if accepts_normal:
                options.append({**kwargs, "normal": True})
            for opts in options:
                for idx, (seq_x, seq_y) in enumerate(PAIRS):
                    key = f"{name}|{sorted(opts.items())}|{idx}"
                    results[key] = _call(func, seq_x, seq_y, opts)
                    results[key + "|rev"] = _call(func, seq_y, seq_x, opts)

    for idx, (seq_x, seq_y) in enumerate(TITLE_PAIRS):
        results[f"alignment.monge_elkan_simil|{idx}"] = _call(
            alignment.monge_elkan_simil, seq_x, seq_y, {}
        )
    for idx, (seq_x, _) in enumerate(PAIRS):
        results[f"compression.lz76_complexity|{idx}"] = repr(
            compression.lz76_complexity(seq_x)
        )
        results[f"order.restrict_to_shared|{idx}"] = repr(
            order.restrict_to_shared(seq_x, PAIRS[idx][1])
        )
        results[f"order.restrict_to_shared|first|{idx}"] = repr(
            order.restrict_to_shared(seq_x, PAIRS[idx][1], repeats="first")
        )

    return results


if __name__ == "__main__":
    SNAPSHOT.parent.mkdir(exist_ok=True)
    SNAPSHOT.write_text(json.dumps(compute(), indent=0, sort_keys=True) + "\n")
    print(f"Wrote {SNAPSHOT}")
