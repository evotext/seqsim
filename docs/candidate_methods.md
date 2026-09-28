# Candidate methods: status

Research notes on measures from the literature considered for `seqsim`, and
what was done with each. Claims about metric properties were checked against
the primary sources and, for the implemented methods, by exhaustive or
property-based tests.

## Implemented (version 0.4.0)

| Method | Function | Notes |
|---|---|---|
| Breakpoint | `order.breakpoint_dissim` | Ordered adjacencies with boundaries; pseudometric with repeated elements |
| IEBP | `order.iebp_estimate` | Formulas of Spencer et al. (2003) for unsigned linear orders under transpositions |
| Kendall tau (top-k) | `order.kendall_tau_dissim` | Fagin et al.'s K^(p) is a near metric, not a metric, for every p |
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
