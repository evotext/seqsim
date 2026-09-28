# Candidate methods for seqsim (research summary)

Compiled by a research subagent; references marked [UNVERIFIED] were not confirmed online.

## Shortlist
1. **Breakpoint / adjacency dissimilarity**: a metric on permutations and a pseudometric in general. O(n). Used for Canterbury Tales tale order (Spencer et al. 2003, doi:10.1023/A:1021818600001; Sankoff & Blanchette 1998, doi:10.1089/cmb.1998.5.555).
2. **Rank distances**: Kendall tau, Spearman footrule, Ulam, Cayley. All are metrics on permutations and need a policy for missing or duplicated items (Diaconis & Graham 1977, doi:10.1111/j.2517-6161.1977.tb01624.x; Aldous & Diaconis 1999).
3. **Top-k / partial-list Kendall and footrule**: handle differing contents. Some variants are metrics; check the exact parameter ranges (Fagin, Kumar & Sivakumar 2003, SIAM J. Discrete Math. 17(1)).
4. **LCS / indel distance** (metric), plus **Bakkelund's normalized LCS metric** (Bakkelund 2009 report).
5. **Metric normalizations of edit distance in [0,1]**: Yujian & Bo GLD (doi:10.1109/TPAMI.2007.1078); Marzal-Vidal, proven a metric for uniform costs by Fisman et al. 2022 (doi:10.4230/LIPIcs.CPM.2022.17).
6. **Block-move measures**: Greedy String Tiling (Wise 1993; similarity) and Tichy block-move cover (doi:10.1145/357401.357404; asymmetric).
7. **Block interchange** (Christie 1996; metric on permutations) and **DCJ / DCJ-indel** (Bergeron et al. 2006; Braga et al. 2011, doi:10.1089/cmb.2011.0118). Avoid transposition distance (NP-hard) and unsigned reversal distance (NP-hard).
8. **Alignment family**: Needleman-Wunsch / Gotoh with a user cost function (doi:10.1016/0022-2836(82)90398-9), Smith-Waterman local similarity (doi:10.1016/0022-2836(81)90087-5), and Monge-Elkan for sequences of sequences.

Runners-up:
- Ukkonen q-gram distance (doi:10.1016/0304-3975(92)90143-4)
- Lempel-Ziv / Otu-Sayood distance (doi:10.1093/bioinformatics/btg295)
- Broder containment (doi:10.1109/SEQUEN.1997.666900) and Tversky index (doi:10.1037/0033-295X.84.4.327)

Out of scope as pairwise sequence measures: DTW and Fréchet (degenerate with categorical elements); Roos & Heikkilä 2009 and the `stemmatology` R package (these work on collation tables).
