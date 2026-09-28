"""
Analysis of traditions: collections of witnesses sharing a reference frame.

A tradition is a set of witnesses, each an ordered sequence of hashable
items (texts, sayings, strophes, chapters), that differ in which items each
witness has and in what order. This subpackage prepares such data for
phylogenetic software:

  * a reference frame (a consensus order of all items) and a label
    (chapter, book) for each item, from the witnesses and their layouts;
  * coverage: where the absence of an item from a witness is evidence,
    and where it is missing data;
  * binary content and adjacency characters;
  * export to PHYLIP and NEXUS files, partitions, and MrBayes lines.

It uses only the Python standard library and does not build trees. A
witness is any sequence of hashable items, a collection is a mapping from
witness names to witnesses, and a character cell is `1`, `0`, or `None`
(missing). All results are deterministic: they do not depend on the order
in which Python iterates sets (`PYTHONHASHSEED`).
"""

from ._frame import (
    consensus_order,
    fill_labels,
    monotone_labels,
    proportional_labels,
    proportional_positions,
    reference_labels,
    unique,
)
from ._coverage import (
    Cell,
    coverage,
    covered_labels,
    density,
    merge_coverage,
    merge_parts,
)
from ._characters import (
    CharacterMatrix,
    adjacencies,
    adjacency_characters,
    character_blocks,
    concat_characters,
    content_characters,
    informative,
    restrict,
)
from ._export import (
    charset_ranges,
    mrbayes_calibrations,
    mrbayes_constraints,
    partition_nexus,
    taxon_name,
    to_nexus,
    to_phylip,
)

__all__ = [
    "Cell",
    "CharacterMatrix",
    "adjacencies",
    "adjacency_characters",
    "character_blocks",
    "charset_ranges",
    "concat_characters",
    "consensus_order",
    "content_characters",
    "coverage",
    "covered_labels",
    "density",
    "fill_labels",
    "informative",
    "merge_coverage",
    "merge_parts",
    "monotone_labels",
    "mrbayes_calibrations",
    "mrbayes_constraints",
    "partition_nexus",
    "proportional_labels",
    "proportional_positions",
    "reference_labels",
    "restrict",
    "taxon_name",
    "to_nexus",
    "to_phylip",
    "unique",
]
