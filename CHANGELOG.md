# Changelog

## Version 0.4.0

This release fixes several correctness bugs and renames most functions so that
their names state their mathematical properties. It is not backwards
compatible; see the migration guide below.

### Added

New measures, most of them suited to comparing the order of texts in
manuscripts with overlapping contents (see the README for their properties):

- `edit`: `indel_dist`, `lcs_dist` (Bakkelund), the metric normalizations
  `levenshtein_gld_dist`, `damerau_gld_dist` and `indel_gld_dist` (Yujian &
  Bo) and `levenshtein_ned_dist` (Marzal & Vidal), `block_move_dissim`
  (Tichy), and `gst_dissim` (Greedy String Tiling).
- `order` (new module): `ulam_dist`, `kendall_tau_dissim` and
  `footrule_dissim` (generalized to different contents following Fagin et
  al.), `cayley_dissim`, `block_interchange_dissim` (Christie),
  `breakpoint_dissim`, and `iebp_estimate` (Wang & Warnow; Spencer et al.).
- `alignment` (new module): `nw_dissim` (global alignment with custom costs
  and affine gaps), `sw_simil` (local alignment), and `monge_elkan_simil`
  (sequences of sequences).
- `token`: `qgram_dissim` (Ukkonen), `tversky_simil`, and `containment`
  (Broder).
- `compression`: `lz76_complexity` and `lz76_dissim` (Otu & Sayood).
- `common.lcs_length`.

New documentation, with a guide to the concepts and to choosing a method,
tutorials, a case study on a real textual tradition, and a page for every
family of measures. `CITATION.cff` added.

### Breaking changes

- Function names now follow a convention: `_dist` for true distances
  (metrics), `_dissim` for dissimilarities (`0.0` for identical sequences,
  but the metric properties are not all guaranteed), and `_simil` for
  similarities.
- All parameters after the two sequences are keyword-only (e.g.
  `bulk_delete_dist(x, y, max_del_len=3)`), and all measures return floats.
- `arith_ncd` was removed; use `lzma_ncd_dissim`. The arithmetic-coding NCD
  used a model that cannot detect repetition, so identical sequences scored
  1.0 or more at any length.
- `damerau` now computes the true (unrestricted) Damerau-Levenshtein
  distance. The previous algorithm, optimal string alignment, is available as
  `osa`.
- Several measures were made symmetric (see "Changed results" below).
- `distance()` returns dissimilarities only, so `birnbaum_simil` is no longer
  in `METHODS`; call `edit.birnbaum_simil()` directly.
- Python 3.10 or later is required. The package has no third-party
  dependencies (`textdistance`, `numpy`, and `tabulate` were removed).

### Migration guide

| 0.3.1 | 0.4.0 | `distance()` key |
|---|---|---|
| `edit.levenshtein_dist` | `edit.levenshtein_dist` | `levenshtein` |
| `edit.levdamerau_dist` | `edit.osa_dissim` (same algorithm) or `edit.damerau_dist` (true Damerau-Levenshtein) | `osa` / `damerau` |
| `edit.bulk_delete_dist` | `edit.bulk_delete_dist` | `bulk_delete` |
| `edit.fragile_ends_simil` | `edit.fragile_ends_dissim` | `fragile_ends` (was `fragile_ends_simil`) |
| `edit.stemmatological_simil` | `edit.stemmatological_dissim` | `stemmatological` (was `stemmatological_simil`) |
| `edit.jaro_dist` | `edit.jaro_dissim` | `jaro` |
| `edit.jaro_winkler_dist` | `edit.jaro_winkler_dissim` | `jaro_winkler` |
| `edit.mmcwpa_dist` | `edit.mmcwpa_dissim` | `mmcwpa` |
| `edit.birnbaum_dist` | `edit.birnbaum_dissim` | `birnbaum` |
| `edit.fast_birnbaum_dist` | removed (identical to `birnbaum_dist`) | |
| `edit.birnbaum_simil` | `edit.birnbaum_simil` | (removed from `METHODS`) |
| `token.jaccard_dist` | `token.jaccard_dissim` | `jaccard` |
| `token.sorensen_dist` | `token.sorensen_dissim` | `sorensen` |
| `token.subseq_jaccard_dist` | `token.subseq_jaccard_dissim` | `subseq_jaccard` |
| `sequence.ratcliff_obershelp` | `sequence.ratcliff_obershelp_dissim` | `ratcliff_obershelp` |
| `compression.entropy_ncd` | `compression.entropy_ncd_dissim` | `entropy_ncd` (was `entropy`) |
| `compression.arith_ncd` | `compression.lzma_ncd_dissim` (different method) | `lzma_ncd` (was `arith_ncd`) |

### Fixed

- `mmcwpa`: only the first subfield was searched, giving wrong and asymmetric
  results (kitten/sitting is now 0.5135, was 0.5385).
- `distance()`: the mean over more than two sequences divided by the number of
  sequences instead of the number of pairs, inflating results for four or more
  sequences.
- `stemmatological`: the largest block deletion was one element shorter than
  `max_del_len`, and `max_del_len=1` disabled interior deletions.
- `bulk_delete` and `stemmatological`: invalid `max_del_len` values (0 or
  negative) crashed or returned negative distances; they now raise
  `ValueError`, as do `frag_start`/`frag_end` values outside [0..100].
- `birnbaum_simil`: identical sequences ignored `normal` and returned an int.
- `ratcliff_obershelp`: elements were joined as strings, so e.g. `[1, 23]` and
  `[12, 3]` were considered identical.
- `subseq_jaccard`: identical sequences with repeated sub-sequences had a
  non-zero score, and results were raised to the power of the sequence length,
  pushing long sequences towards zero.
- Empty sequences raised `ZeroDivisionError` in many methods. Two empty
  sequences now score `0.0`, and an empty sequence compared with a non-empty
  one scores the maximum (`1.0` in range [0..1]).
- `difflib`'s automatic junk heuristic no longer changes results for
  sequences of 200 or more elements.
- n-gram padding uses a dedicated `PAD` sentinel instead of the string
  `"$$$"`, which could collide with sequence elements.
- The source distribution could not be installed (`requirements.txt` was not
  shipped).

### Changed results

- `bulk_delete`: blocks can be inserted as well as deleted, making it
  symmetric and a true metric.
- `fragile_ends` and `stemmatological`: the discount applies to insertions
  and deletions at the fragile ends of both sequences (symmetric); a block is
  discounted only if it lies entirely within a fragile region; the end region
  no longer includes one extra position; `fragile_ends` is normalized by the
  longest length, like the other edit measures.
- `jaro`, `jaro_winkler`, `ratcliff_obershelp`, `mmcwpa`, `birnbaum_simil`:
  the greedy matching is computed in both argument orders and the most
  similar result is used.
- `birnbaum`: normalized by the longest sequence (it was the shortest, so any
  sequence contained in another scored `0.0`).
- `subseq_jaccard`: multiset Jaccard per sub-sequence length, averaged with
  length weights (no exponent).
- `mmcwpa` and `subseq_jaccard` are much faster (about 300x and 5x for 400
  elements).

## Version 0.3.1

- Fixed bug due to typo in one of the methods
- Selected one Birnbaum implementation

## Version 0.3

- Improvements to code quality, documentation, and references
- Added new methods and scaffolding for future expansions

## Version 0.2

- First release for new roadmap supporting sequences of any hashable Python
  datatype, importing code from other projects (mostly from `titivillus`)
