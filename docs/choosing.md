# Choosing a method

There is no single best measure: each one models a different idea of what
makes two sequences "different". This page starts from common research
questions and points to the measures that answer them. The
[method pages](methods/index.md) give the details of each one.

## Quick guide

| Question | Elements | Recommended | Alternatives |
|---|---|---|---|
| How many variant readings separate two witnesses? | words | `edit.levenshtein_dist` | `edit.indel_dist`, `alignment.nw_dissim` |
| The same, as a proportion, for building a stemma or tree | words | `edit.levenshtein_gld_dist` | `edit.levenshtein_ned_dist`, `edit.lcs_dist` |
| Are transpositions of words common in the tradition? | words | `edit.damerau_dist` | `edit.damerau_gld_dist` |
| How different are two spellings? | characters | `edit.levenshtein_dist` (normalized) | `edit.jaro_winkler_dissim`, `alignment.nw_dissim` with custom costs |
| Which texts do two manuscripts share? | texts | `token.jaccard_dissim` | `token.sorensen_dissim` |
| Is one manuscript an excerpt of another? | texts or words | `token.containment` | `token.tversky_simil` |
| How different is the order of the shared texts? | texts | `order.ulam_dist` | `order.kendall_tau_dissim`, `order.breakpoint_dissim` |
| How many blocks of texts (quires, gatherings) were moved? | texts | `order.block_interchange_dissim` | `edit.block_move_dissim`, `edit.gst_dissim` |
| How many rearrangements happened, correcting for multiple changes? | texts | `order.iebp_estimate` | |
| Were texts lost at the beginning or end (damaged manuscripts)? | texts | `edit.stemmatological_dissim` | `edit.fragile_ends_dissim`, `edit.bulk_delete_dist` |
| Were whole groups of texts added or lost? | texts | `edit.bulk_delete_dist` | `edit.stemmatological_dissim` |
| Do two long texts share passages, in any order? | words | `edit.gst_dissim` | `alignment.sw_simil`, `token.qgram_dissim` |
| Which titles or incipits refer to the same text? | characters | `alignment.monge_elkan_simil` | `edit.jaro_winkler_dissim` |
| A quick, alignment-free comparison of long texts | words | `token.qgram_dissim` | `compression.lzma_ncd_dissim`, `compression.lz76_dissim` |

## Comparing witnesses of a text

When the sequences are the words of different witnesses of the same text, the
natural measures are **edit distances**, which count the variant readings
(substitutions, additions, and omissions) needed to turn one witness into the
other.

- `edit.levenshtein_dist` counts substitutions, insertions and deletions of
  single words. It is the standard choice.
- `edit.indel_dist` counts only additions and omissions, so that a
  substitution counts as two changes; `edit.lcs_dist` is its normalized
  counterpart, based on the longest common subsequence.
- `edit.damerau_dist` also counts the transposition of two adjacent words as
  a single change, which is common in some traditions.
- `alignment.nw_dissim` lets you define the cost of each substitution, for
  example to make orthographic variants cheaper than substantive ones, and to
  make a long omission (an eye-skip, or *saut du même au même*) cost less than
  the same number of isolated omissions, with affine gaps.

For a matrix that will be used to build a tree or network, use a normalized
measure that is also a true distance, such as `edit.levenshtein_gld_dist`, so
that witnesses of different lengths are comparable (see
[Normalization](concepts.md#normalization)).

The choice of elements matters as much as the measure. Comparing raw words
counts spelling variants as differences; normalizing the spelling (or
comparing lemmata) before the comparison focuses on substantive variants.

## Comparing the contents of manuscripts

When the sequences are the lists of texts copied in manuscripts (for example,
miscellanies, collections of saints' lives, sermons, or songs), two questions
are usually separated:

**Which texts are shared?** Set-based measures, such as
`token.jaccard_dissim`, ignore the order and measure the overlap of the
contents. `token.containment` tells how much of one manuscript is found in
another, which helps detecting excerpts and derived collections.

**How different is the order of the shared texts?** The measures of the
`order` module come from the comparison of rankings and from the study of
genome rearrangements, where the same questions arise for genes on
chromosomes:

- `order.ulam_dist` counts the texts that must be moved, one at a time
  (texts present in only one manuscript count as an insertion or deletion).
- `order.kendall_tau_dissim` counts the pairs of texts in a different
  relative order, and `order.footrule_dissim` sums how far each text moved.
- `order.block_interchange_dissim` counts exchanges of two blocks of texts,
  modelling quires or gatherings bound in a different order.
- `order.breakpoint_dissim` counts the pairs of consecutive texts that are
  no longer consecutive, the measure used in the study of the order of the
  *Canterbury Tales* by Spencer et al. (2003), and `order.iebp_estimate`
  corrects it into an estimate of the number of rearrangements.

**What was lost or added, and where?** Physical damage often removes the
first or last texts of a manuscript, and whole gatherings can be lost or
added. The "stemmatological" dissimilarity (`edit.stemmatological_dissim`)
was designed for this: it counts the loss of a block of consecutive texts as
a single event, and discounts losses at the beginning and end of the
manuscript.

The [manuscript contents tutorial](tutorials/contents.md) walks through these
measures with an example.

## Identifying texts

Before comparing contents, texts must be identified across manuscripts, often
from titles, rubrics, or incipits with variable spelling and wording.
`alignment.monge_elkan_simil` compares two lists of titles word by word,
tolerating differences in word order and spelling; `edit.jaro_winkler_dissim`
and the normalized Levenshtein distance compare single titles. See the
[tutorial on identifying texts](tutorials/identifying.md).

## When in doubt

- Prefer simple, well-understood measures (Levenshtein, Jaccard, Ulam), and
  use the others to answer specific questions.
- Prefer a `_dist` measure when the results will be used to infer a tree or a
  network.
- Compare the results of two or three measures: robust conclusions should not
  depend on the choice of one of them.
