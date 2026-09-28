# Methods

The measures are organized in modules by family. Each method page explains
what every measure models, how it is computed, its mathematical properties,
its parameters, and when to use it, with references to the original
publications.

- [Edit measures](edit.md) (`seqsim.edit`): edit distances, their
  normalizations, block operations, and measures developed for manuscript
  traditions.
- [Order measures](order.md) (`seqsim.order`): the order of shared items,
  from the comparison of rankings and genome rearrangements.
- [Alignment](alignment.md) (`seqsim.alignment`): global and local alignment
  with custom costs, and the comparison of lists of titles.
- [Token measures](token.md) (`seqsim.token`): shared elements and
  sub-sequences, regardless of their position.
- [Sequence matching](sequence.md) (`seqsim.sequence`): Ratcliff-Obershelp.
- [Compression](compression.md) (`seqsim.compression`): measures based on
  compression and on the complexity of sequences.
- [Traditions](tradition.md) (`seqsim.tradition`): reference frames,
  coverage, and content and adjacency characters for whole collections of
  witnesses, with export for phylogenetic software.

## Overview

The table lists the measures available through `seqsim.distance()`, with
their properties (see [Concepts](../concepts.md)): whether only identical
sequences score zero (identity), and whether the triangle inequality holds.
All of them are symmetric and accept `normal=True`.

| `distance()` key | Function | Identity | Triangle inequality | Raw range |
|---|---|---|---|---|
| `levenshtein` | [`edit.levenshtein_dist`](edit.md#levenshtein) | yes | yes | 0 to max length |
| `damerau` | [`edit.damerau_dist`](edit.md#damerau-levenshtein) | yes | yes | 0 to max length |
| `indel` | [`edit.indel_dist`](edit.md#indel-and-lcs) | yes | yes | 0 to sum of lengths |
| `lcs` | [`edit.lcs_dist`](edit.md#indel-and-lcs) | yes | yes | 0 to 1 |
| `levenshtein_gld` | [`edit.levenshtein_gld_dist`](edit.md#normalized-edit-distances) | yes | yes | 0 to 1 |
| `damerau_gld` | [`edit.damerau_gld_dist`](edit.md#normalized-edit-distances) | yes | yes | 0 to 1 |
| `indel_gld` | [`edit.indel_gld_dist`](edit.md#normalized-edit-distances) | yes | yes | 0 to 1 |
| `levenshtein_ned` | [`edit.levenshtein_ned_dist`](edit.md#normalized-edit-distances) | yes | yes | 0 to 1 |
| `bulk_delete` | [`edit.bulk_delete_dist`](edit.md#bulk-delete) | yes | yes | 0 to max length |
| `ulam` | [`order.ulam_dist`](order.md#ulam) | yes | yes | 0 to sum of lengths |
| `osa` | [`edit.osa_dissim`](edit.md#optimal-string-alignment) | yes | no | 0 to max length |
| `fragile_ends` | [`edit.fragile_ends_dissim`](edit.md#fragile-ends) | yes | no | 0 to max length |
| `stemmatological` | [`edit.stemmatological_dissim`](edit.md#stemmatological) | yes | no | 0 to max length |
| `block_move` | [`edit.block_move_dissim`](edit.md#block-moves) | yes | no | 0 to max length + 1 |
| `gst` | [`edit.gst_dissim`](edit.md#greedy-string-tiling) | no | no | 0 to 1 |
| `jaro` | [`edit.jaro_dissim`](edit.md#jaro-and-jaro-winkler) | yes | no | 0 to 1 |
| `jaro_winkler` | [`edit.jaro_winkler_dissim`](edit.md#jaro-and-jaro-winkler) | yes | no | 0 to 1 |
| `mmcwpa` | [`edit.mmcwpa_dissim`](edit.md#mmcwpa) | yes | no | 0 to 1 |
| `birnbaum` | [`edit.birnbaum_dissim`](edit.md#birnbaum) | yes | no | 0 to 1 |
| `kendall_tau` | [`order.kendall_tau_dissim`](order.md#kendall-tau) | yes | no | 0 to number of pairs |
| `footrule` | [`order.footrule_dissim`](order.md#spearman-footrule) | yes | with a fixed `ell` | 0 upwards |
| `cayley` | [`order.cayley_dissim`](order.md#cayley) | yes | no | 0 to sum of lengths |
| `block_interchange` | [`order.block_interchange_dissim`](order.md#block-interchange) | yes | not proven | 0 to sum of lengths |
| `breakpoint` | [`order.breakpoint_dissim`](order.md#breakpoints) | no | yes | 0 upwards |
| `nw` | [`alignment.nw_dissim`](alignment.md#global-alignment) | yes | depends on costs | 0 upwards |
| `jaccard` | [`token.jaccard_dissim`](token.md#jaccard) | no | yes | 0 to 1 |
| `sorensen` | [`token.sorensen_dissim`](sorensen-dice) | no | no | 0 to 1 |
| `subseq_jaccard` | [`token.subseq_jaccard_dissim`](token.md#sub-sequence-jaccard) | yes | not proven | 0 to 1 |
| `qgram` | [`token.qgram_dissim`](token.md#q-grams) | no | yes | 0 upwards |
| `ratcliff_obershelp` | [`sequence.ratcliff_obershelp_dissim`](sequence.md) | yes | no | 0 to 1 |
| `entropy_ncd` | [`compression.entropy_ncd_dissim`](compression.md#entropy-ncd) | no | no | 0 to 1 |
| `lzma_ncd` | [`compression.lzma_ncd_dissim`](compression.md#lzma-ncd) | no | not guaranteed | about 0 to 1 |
| `lz76` | [`compression.lz76_dissim`](compression.md#lempel-ziv) | no | no | about 0 to 1 |

Other functions:

| Function | Kind |
|---|---|
| [`edit.birnbaum_simil`](edit.md#birnbaum) | similarity |
| [`alignment.sw_simil`](alignment.md#local-alignment) | similarity (local alignment) |
| [`alignment.monge_elkan_simil`](alignment.md#monge-elkan) | similarity of sequences of sequences |
| [`token.tversky_simil`](token.md#tversky) | similarity, directional unless `alpha == beta` |
| [`token.containment`](token.md#containment) | directional proportion |
| [`order.kendall_tau_simil`](order.md#kendall-correlation-of-the-shared-items) | correlation of the orders of the shared items, in range [-1..1] |
| [`order.breakpoint_simil`](order.md#breakpoints) | share of adjacencies preserved |
| [`order.restrict_to_shared`](order.md#restricting-to-the-shared-items) | reduces two sequences to their shared items |
| [`order.iebp_estimate`](order.md#iebp) | estimate of the number of rearrangements |
| [`compression.lz76_complexity`](compression.md#lempel-ziv) | complexity of a single sequence |
| `common.lcs_length` | length of the longest common subsequence |
| `common.equivalent_string` | maps two sequences to equivalent strings |
