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

```{include} _measures.md
```

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
