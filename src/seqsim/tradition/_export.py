"""
Export of character matrices for phylogenetic software.

The functions write relaxed PHYLIP (for IQ-TREE and RAxML), NEXUS (for
MrBayes and PAUP*), the NEXUS `sets` blocks that define partitions, and the
MrBayes lines for tip dates and constrained, calibrated groups. They encode
no modelling choice: models and priors are left to the user.
"""

# Import Python standard libraries
from typing import (
    Callable,
    Dict,
    Hashable,
    List,
    Mapping,
    Optional,
    Sequence,
    Tuple,
)

# Import local modules
from ._characters import CharacterMatrix

_SYMBOL = {1: "1", 0: "0"}

TAXON_REPLACEMENTS: Tuple[Tuple[str, str], ...] = ((" ", "_"), ("+", "-"))


def taxon_name(
    witness_id: str,
    replacements: Sequence[Tuple[str, str]] = TAXON_REPLACEMENTS,
) -> str:
    """
    Makes a witness name safe for PHYLIP, NEXUS, Newick, and IQ-TREE.

    By default, spaces become underscores and the `+` joining merged parts
    becomes `-`. Other replacements, applied in order, can be given for
    names with other problematic characters (such as parentheses, colons,
    or commas).

    Example
    ********

    .. code-block:: python

        >>> seqsim.tradition.taxon_name("Paris gr. 1596 A+C")
        'Paris_gr._1596_A-C'
        >>> seqsim.tradition.taxon_name("Vat. gr. 1 (A)", [(" ", "_"), ("(", ""), (")", "")])
        'Vat._gr._1_A'

    :param witness_id: The name of the witness.
    :param replacements: A sequence of `(old, new)` string replacements.
    :return: The taxon name.
    """

    name = witness_id
    for old, new in replacements:
        name = name.replace(old, new)

    return name


def _states(row: Sequence) -> str:
    return "".join(_SYMBOL.get(cell, "?") for cell in row)


def to_phylip(
    matrix: CharacterMatrix, rename: Callable[[str], str] = taxon_name
) -> str:
    """
    Writes a character matrix in relaxed PHYLIP format, with `?` for missing.

    The first line gives the number of taxa and of characters; each taxon
    follows on its own line, with its name (see `taxon_name`), a space, and
    its states.

    Example
    ********

    .. code-block:: python

        >>> matrix = seqsim.tradition.CharacterMatrix(
        ...     ("W 1", "W 2"), ("a", "b", "c"), ((1, 0, None), (0, 1, 1))
        ... )
        >>> print(seqsim.tradition.to_phylip(matrix), end="")
        2 3
        W_1 10?
        W_2 011

    :param matrix: The character matrix.
    :param rename: The function giving the name of each taxon. Defaults to
        `taxon_name`.
    :return: The PHYLIP text, ending with a newline.
    """

    lines = [f"{matrix.shape[0]} {matrix.shape[1]}"]
    for taxon, row in zip(matrix.taxa, matrix.cells):
        lines.append(f"{rename(str(taxon))} " + _states(row))

    return "\n".join(lines) + "\n"


def to_nexus(
    matrix: CharacterMatrix,
    charsets: Optional[Mapping[str, Tuple[int, int]]] = None,
    rename: Callable[[str], str] = taxon_name,
) -> str:
    """
    Writes a character matrix in NEXUS format, for restriction-site models.

    The matrix is written as a `data` block with `datatype=restriction` (the
    binary type of MrBayes), `?` for missing, and the names aligned, followed
    by an optional `sets` block with a `charset` for each partition.

    Example
    ********

    .. code-block:: python

        >>> matrix = seqsim.tradition.CharacterMatrix(
        ...     ("W1", "W2"), ("a", "b"), ((1, None), (0, 1))
        ... )
        >>> print(seqsim.tradition.to_nexus(matrix, {"content": (1, 2)}), end="")
        #NEXUS
        begin data;
          dimensions ntax=2 nchar=2;
          format datatype=restriction missing=? gap=-;
          matrix
            W1  1?
            W2  01
          ;
        end;
        begin sets;
          charset content = 1-2;
        end;

    :param matrix: The character matrix.
    :param charsets: An optional mapping from partition names to the first
        and last character (counting from 1) of each partition.
    :param rename: The function giving the name of each taxon. Defaults to
        `taxon_name`.
    :return: The NEXUS text, ending with a newline.
    """

    names = [rename(str(taxon)) for taxon in matrix.taxa]
    width = max((len(name) for name in names), default=0) + 2
    lines = [
        "#NEXUS",
        "begin data;",
        f"  dimensions ntax={matrix.shape[0]} nchar={matrix.shape[1]};",
        "  format datatype=restriction missing=? gap=-;",
        "  matrix",
    ]
    for name, row in zip(names, matrix.cells):
        lines.append(f"    {name.ljust(width)}" + _states(row))
    lines += ["  ;", "end;"]
    if charsets:
        lines += ["begin sets;"]
        lines += [f"  charset {name} = {a}-{b};" for name, (a, b) in charsets.items()]
        lines += ["end;"]

    return "\n".join(lines) + "\n"


def charset_ranges(sizes: Mapping[str, int]) -> Dict[str, Tuple[int, int]]:
    """
    Returns the ranges of consecutive partitions, from their sizes.

    Example
    ********

    .. code-block:: python

        >>> seqsim.tradition.charset_ranges({"content": 120, "adjacency": 80})
        {'content': (1, 120), 'adjacency': (121, 200)}

    :param sizes: A mapping from partition names to their number of
        characters, in the order of the partitions in the matrix.
    :return: A mapping from partition names to their first and last
        character, counting from 1.
    """

    ranges = {}
    start = 1
    for name, size in sizes.items():
        ranges[name] = (start, start + size - 1)
        start += size

    return ranges


def partition_nexus(charsets: Mapping[str, Tuple[int, int]]) -> str:
    """
    Writes a NEXUS `sets` block defining partitions, as read by IQ-TREE.

    Example
    ********

    .. code-block:: python

        >>> charsets = seqsim.tradition.charset_ranges({"content": 120, "adjacency": 80})
        >>> print(seqsim.tradition.partition_nexus(charsets), end="")
        #nexus
        begin sets;
          charset content = 1-120;
          charset adjacency = 121-200;
        end;

    :param charsets: A mapping from partition names to their first and last
        character, counting from 1 (see `charset_ranges`).
    :return: The NEXUS text, ending with a newline.
    """

    lines = ["#nexus", "begin sets;"]
    lines += [f"  charset {name} = {a}-{b};" for name, (a, b) in charsets.items()]
    lines += ["end;"]

    return "\n".join(lines) + "\n"


def mrbayes_calibrations(
    tip_ages: Mapping[str, Tuple[Hashable, Hashable]],
) -> List[str]:
    """
    Writes the MrBayes `calibrate` lines for the ages of the tips.

    A tip whose minimum and maximum ages are equal gets a fixed age, and the
    others a uniform prior between the two. Ages are written as given, so
    pass them with the type (integer or float) wanted in the output.

    Example
    ********

    .. code-block:: python

        >>> ages = {"W1": (500, 500), "W2": (300, 400)}
        >>> for line in seqsim.tradition.mrbayes_calibrations(ages):
        ...     print(line)
          calibrate W1 = fixed(500);
          calibrate W2 = uniform(300,400);

    :param tip_ages: A mapping from taxon names to their minimum and maximum
        ages.
    :return: A list of lines, indented by two spaces, for a `mrbayes` block.
    """

    lines = []
    for taxon, (age_min, age_max) in tip_ages.items():
        if age_min == age_max:
            lines.append(f"  calibrate {taxon} = fixed({age_min});")
        else:
            lines.append(f"  calibrate {taxon} = uniform({age_min},{age_max});")

    return lines


def mrbayes_constraints(
    groups: Mapping[str, Tuple[Sequence[str], Optional[Tuple[Hashable, Hashable]]]],
) -> List[str]:
    """
    Writes the MrBayes lines constraining (and dating) groups of taxa.

    Each group is given by its members and, optionally, the youngest and
    oldest age of its common ancestor. The lines define a `constraint` for
    each group, a uniform `calibrate` for each dated group, and a
    `prset topologypr` enforcing all constraints. The node age prior itself
    (`prset nodeagepr=calibrated`) and the other priors are left to the
    user.

    Example
    ********

    .. code-block:: python

        >>> groups = {"west": (["W1", "W2"], (600, 900)), "east": (["W3", "W4"], None)}
        >>> for line in seqsim.tradition.mrbayes_constraints(groups):
        ...     print(line)
          constraint west = W1 W2;
          constraint east = W3 W4;
          calibrate west = uniform(600,900);
          prset topologypr=constraints(west,east);

    :param groups: A mapping from group names to their members and to the
        youngest and oldest age of their ancestor (or `None`).
    :return: A list of lines, indented by two spaces, for a `mrbayes` block;
        empty if there are no groups.
    """

    if not groups:
        return []

    lines = []
    for name, (members, _) in groups.items():
        lines.append(f"  constraint {name} = {' '.join(members)};")
    for name, (_, ages) in groups.items():
        if ages is not None:
            young, old = ages
            lines.append(f"  calibrate {name} = uniform({young},{old});")
    lines.append(f"  prset topologypr=constraints({','.join(groups)});")

    return lines
