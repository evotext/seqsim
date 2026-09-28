# Candidate methods: status

Research notes on measures from the literature considered for `seqsim`, and
what was done with each. Claims about metric properties were checked against
the primary sources and, for the implemented methods, by exhaustive or
property-based tests.

## Implemented (version 0.4.0)

| Method | Function | Notes |
|---|---|---|
| Breakpoint | `order.breakpoint_dissim`, `order.breakpoint_simil` | Ordered adjacencies, with boundaries by default; pseudometric with repeated elements |
| IEBP | `order.iebp_estimate` | Formulas of Spencer et al. (2003) for unsigned linear orders under transpositions |
| Kendall tau (top-k) | `order.kendall_tau_dissim` | Fagin et al.'s K^(p) is a near metric, not a metric, for every p |
| Kendall correlation | `order.kendall_tau_simil` | Over the shared items only; equals tau-a and tau-b without ties |
| Spearman footrule (top-k) | `order.footrule_dissim` | Metric with a fixed `ell` (Fagin et al., Prop. 3.8) |
| Ulam | `order.ulam_dist` | Generalized to moves, insertions and deletions: `|x| + |y| - M - LCS` |
| Cayley | `order.cayley_dissim` | Not a metric with repeated elements |
| Block interchange | `order.block_interchange_dissim` | Christie (1996), `(n + 1 - c) / 2` |
| Indel / LCS | `edit.indel_dist`, `edit.lcs_dist` | Bakkelund (2009), Theorem 3.1 |
| Normalized edit distances | `edit.*_gld_dist`, `edit.levenshtein_ned_dist` | Yujian & Bo (2007); Marzal & Vidal (1993), metric for unit costs (Fisman et al. 2022, Thm. 27) |
| Tichy block moves | `edit.block_move_dissim` | Greedy minimal cover (Tichy 1984) |
| Greedy String Tiling | `edit.gst_dissim` | Wise (1993), as described in JPlag (Prechelt et al. 2002) |
| Needleman-Wunsch / Gotoh | `alignment.nw_dissim` | |
| Smith-Waterman | `alignment.sw_simil` | |
| Monge-Elkan | `alignment.monge_elkan_simil` | |
| q-gram distance | `token.qgram_dissim` | Ukkonen (1992); padded by default |
| Tversky, containment | `token.tversky_simil`, `token.containment` | Directional |
| Lempel-Ziv | `compression.lz76_dissim` | Otu & Sayood (2003), measure d* |

## Measures of the Apophthegmata project (version 0.4.0)

The research project on the *Apophthegmata Patrum*, for which `seqsim` is
the companion library, used three pairwise measures of its own. Each was
compared with its `seqsim` counterpart on all pairs of its 48 witnesses, and
the following decisions were taken:

| Project measure | Decision |
|---|---|
| `containment(a, b)`, the share of the distinct items of `a` found in `b` | Kept `token.containment` (`size=1`), which counts repeats as a multiset. Identical results on sequences reduced to first occurrences (`tradition.unique`), which is how the project calls it; up to 0.05 apart otherwise. |
| `order_tau(a, b)`, Kendall's tau of the shared items (SciPy) | Kept the discordant-pair count of `kendall_tau_dissim`, and added `order.kendall_tau_simil`, the correlation `1 - 4d / (n(n-1))` over the shared items, with `nan` for fewer than two. Identical to SciPy up to 3e-16. A minimum number of shared items is left to the caller, not added as a parameter. |
| `adjacency_agreement(a, b)`, the share of adjacencies of the shared items preserved, without boundaries | Kept the breakpoint of `breakpoint_dissim`, and added a `boundaries` option: without boundaries, a sequence of `n` items has `n - 1` adjacencies, which suits fragments, starting and ending anywhere, and matches the adjacency characters of `tradition`. `order.breakpoint_simil` computes the share directly, so that it is bit-for-bit identical to the project's value (`1 - breakpoint_dissim(normal=True)` can differ in the last bit). |

The restriction to shared items is public as `order.restrict_to_shared`.
Repeated items are matched by occurrence by default, as in all measures of
the `order` module; `repeats="first"` first reduces each sequence to the
first occurrence of each item, as the project does.

The IEBP estimator also accepts `boundaries=False`, summing only the terms of
the formulas of Spencer et al. (2003) for the interior positions; it remains
unbiased for few transpositions in simulations. Without boundaries, the
breakpoint dissimilarity keeps its properties: it is symmetric and satisfies
the triangle inequality (half the size of the symmetric difference of two
multisets), and for permutations of the same (two or more) items only
identical orders score zero; `"a"` and `"b"`, with no adjacencies, score
zero.

The comparison also exposed a problem in the normalization of
`kendall_tau_dissim`, which divided by the number of pairs of items
including the end boundary, `C(n + 1, 2)` for permutations of `n` items. The
pairs of the boundary with an item present in both sequences always cost
zero, so two reversed permutations scored at most `(n - 1) / (n + 1)` (0.6
for four items), and 1.0 was reached only for degenerate pairs such as `""`
and `"a"`. It now divides by its maximum over all orders of the same
contents: all pairs, minus those with the boundary and a shared item, with
the pairs of two items missing from the same sequence counted as `p`. For
permutations this is `C(n, 2)`, the classic normalization. The maximum is
attained (items present in only one sequence first, shared items reversed),
which was checked exhaustively for all pairs of contents of up to four
items; with repeated items, it is an upper bound.

## Not implemented

- **DCJ / DCJ-indel** (Bergeron et al. 2006; Braga et al. 2011): requires an
  orientation of each item, which texts in manuscripts do not have.
- **Transposition distance**: NP-hard (Bulteau, Fertin & Rusu 2012).
- **Unsigned reversal distance**: NP-hard; reversals are also rare for text
  order.
- **Edit distance with block moves** (Cormode & Muthukrishnan; Shapira &
  Storer): NP-hard exactly; `block_move_dissim` and `gst_dissim` cover the
  use case.
- **Dynamic time warping, discrete Fréchet**: degenerate for categorical
  elements.
- **Collation-based stemmatological measures** (e.g., Roos & Heikkilä 2009;
  the `stemmatology` R package): these operate on tables of variant
  locations, not on pairs of sequences.
