"""
The measure contract: declaration, conventions, and registry of all measures.

Every measure of the library is a function decorated with `measure`, which
declares its metadata (its kind, its claims, how it is normalized) in one
place and applies the conventions shared by all measures:

  * parameter checks, before anything else;
  * the empty rule: two empty sequences are identical, and an empty
    sequence compared with a non-empty one has the largest value;
  * symmetrization, computing the measure in both orders of the arguments;
  * normalization, dividing the raw value by a bound when `normal=True`;
  * a float result.

The decorated body only computes the raw value, for one order of the
arguments. The registry built by the decorator is the single source for
`seqsim.METHODS`, `seqsim.measures()`, the property tests, and the tables of
the documentation. See `CONTEXT.md` for the vocabulary.
"""

# Import Python standard libraries
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, NamedTuple, Optional, Tuple
import functools
import inspect
import textwrap

# Import local modules
from .common import empty_dissim

# Kinds of measures, by the suffix of their name (or their role)
KINDS = ("dist", "dissim", "simil", "directional", "estimator")

# Claims about the triangle inequality
TRIANGLE = ("yes", "no", "unproven", "conditional")

# How raw values are normalized when `normal=True`
BOUNDS = {
    "max_len": "the length of the longest sequence",
    "sum_len": "the sum of the lengths of both sequences",
    "unit": None,  # the raw value is already in range [0..1]
    "clip": "clipping the raw value to range [0..1]",
    "scored": None,  # the body returns a `Scored` value with its bound
    "custom": None,  # the body handles `normal` itself
}

# How the values of both orders of the arguments are combined
SYMMETRIZE = {
    "min": min,
    "max": max,
    "mean": lambda value_x, value_y: (value_x + value_y) / 2.0,
}


class Scored(NamedTuple):
    """
    A raw value returned by a measure body together with its bound.

    Used when the bound depends on quantities computed by the algorithm, so
    that it is not computed twice. A bound of zero normalizes to `0.0`.
    """

    value: float
    bound: float


@dataclass(frozen=True)
class MeasureInfo:
    """
    The declaration of a measure.

    :ivar func: The public function.
    :ivar name: The qualified name of the function (e.g. `"edit.levenshtein_dist"`).
    :ivar key: The name of the method in `seqsim.distance()` and `METHODS`,
        or `None` if the measure is not available there.
    :ivar kind: One of `"dist"`, `"dissim"`, `"simil"`, `"directional"`, or
        `"estimator"`.
    :ivar identity: Whether a value of zero (for distances and
        dissimilarities) implies identical sequences: `"yes"`, `"no"`, or
        `"unproven"`.
    :ivar identity_example: For `identity="no"`, two different sequences
        with a value of zero.
    :ivar identical_zero: Whether identical sequences always have a value
        of zero (for similarities, the largest value).
    :ivar symmetric: Whether the measure is symmetric.
    :ivar triangle: The triangle inequality claim for raw values: `"yes"`,
        `"no"`, `"unproven"`, or `"conditional"`.
    :ivar triangle_example: For `triangle="no"`, three sequences `(x, y, z)`
        with `d(x, z) > d(x, y) + d(y, z)`.
    :ivar condition: For `triangle="conditional"`, when the inequality holds.
    :ivar bound: How the measure is normalized (see `BOUNDS`).
    :ivar raw_range: A description of the range of raw values.
    :ivar empty: The empty rule: `"max"` if an empty sequence compared with
        a non-empty one has the largest value (1.0 for dissimilarities, 0.0
        for similarities), or `"natural"` if the algorithm defines it.
    :ivar normal_effect: Whether `normal=True` changes the result.
    :ivar accepts_normal: Whether the function accepts `normal`.
    """

    func: Callable
    name: str
    key: Optional[str]
    kind: str
    identity: Optional[str]
    identity_example: Optional[Tuple[Any, Any]]
    identical_zero: bool
    symmetric: bool
    triangle: Optional[str]
    triangle_example: Optional[Tuple[Any, Any, Any]]
    condition: Optional[str]
    bound: str
    raw_range: str
    empty: str
    normal_effect: bool
    accepts_normal: bool


_REGISTRY: List[MeasureInfo] = []


def _describe_normal(bound: str, normal_doc: Optional[str]) -> str:
    """
    Returns the documentation of the `normal` parameter for a bound.
    """

    if normal_doc:
        return normal_doc
    if bound in ("unit",):
        return "Ignored, as results are always in range [0..1]."
    if bound == "clip":
        return "Whether to clip the result to the range [0..1]."

    return (
        "Whether to normalize the result in range [0..1] by dividing it by "
        f"{BOUNDS[bound]}."
    )


def _document_normal(doc: str, text: str) -> str:
    """
    Inserts the documentation of `normal` before the `:return:` of a docstring.
    """

    marker = "    :return:"
    line = "\n".join(
        textwrap.wrap(
            f":param normal: {text}",
            width=84,
            initial_indent="    ",
            subsequent_indent="        ",
        )
    )
    line += "\n"
    if marker in doc:
        return doc.replace(marker, line + marker, 1)

    return doc + "\n" + line


def measure(
    *,
    key: Optional[str] = None,
    kind: str,
    identity: Optional[str] = "yes",
    identity_example: Optional[Tuple[Any, Any]] = None,
    identical_zero: bool = True,
    symmetric: bool = True,
    triangle: Optional[str] = None,
    triangle_example: Optional[Tuple[Any, Any, Any]] = None,
    condition: Optional[str] = None,
    bound: str = "unit",
    raw_range: Optional[str] = None,
    empty: str = "natural",
    symmetrize: Optional[str] = None,
    directional_option: bool = False,
    check: Optional[Callable[..., None]] = None,
    normal: bool = True,
    normal_doc: Optional[str] = None,
) -> Callable[[Callable], Callable]:
    """
    Declares a measure and applies the shared conventions.

    The decorated body takes the two sequences and the keyword arguments of
    the measure (but not `normal`, unless `bound="custom"`, nor
    `directional`), and returns the raw value for that order of the
    arguments, or a `Scored` value when `bound="scored"`.

    :param key: The name of the method in `seqsim.distance()`, if any.
    :param kind: The kind of measure (see `KINDS`).
    :param identity: The identity of indiscernibles claim (`"yes"`,
        `"no"`, or `"unproven"`), for distances and dissimilarities.
    :param identity_example: Two different sequences with a value of zero,
        required when `identity="no"`.
    :param identical_zero: Whether identical sequences have a value of zero.
    :param symmetric: Whether the measure is symmetric.
    :param triangle: The triangle inequality claim (see `TRIANGLE`), for
        distances and dissimilarities.
    :param triangle_example: A counterexample, required when
        `triangle="no"`.
    :param condition: When the inequality holds, required when
        `triangle="conditional"`.
    :param bound: How raw values are normalized (see `BOUNDS`).
    :param raw_range: A description of the range of raw values; defaults to
        one derived from `bound`.
    :param empty: `"max"` to apply the empty rule before calling the body
        (0.0 for two empty sequences and 1.0 otherwise, for dissimilarities;
        the opposite for similarities), or `"natural"` to leave empty
        sequences to the body.
    :param symmetrize: How to combine the values of both orders of the
        arguments (see `SYMMETRIZE`), or `None` if the body is symmetric.
    :param directional_option: Whether to offer a `directional` parameter,
        returning the value for the given order of the arguments only.
    :param check: A function called with the keyword arguments of the
        measure, to validate them, before anything else.
    :param normal: Whether the measure accepts `normal`.
    :param normal_doc: The documentation of `normal`, when the default one
        derived from `bound` does not apply.
    :return: The decorator.
    """

    if kind not in KINDS:
        raise ValueError(f"Unknown kind {kind!r}.")
    if bound not in BOUNDS:
        raise ValueError(f"Unknown bound {bound!r}.")
    if triangle is not None and triangle not in TRIANGLE:
        raise ValueError(f"Unknown triangle claim {triangle!r}.")
    if triangle == "no" and triangle_example is None:
        raise ValueError("A triangle claim of 'no' needs a counterexample.")
    if triangle == "conditional" and not condition:
        raise ValueError("A conditional triangle claim needs a condition.")
    if identity == "no" and identity_example is None:
        raise ValueError("An identity claim of 'no' needs an example.")
    if empty not in ("max", "natural"):
        raise ValueError(f"Unknown empty rule {empty!r}.")
    if symmetrize is not None and symmetrize not in SYMMETRIZE:
        raise ValueError(f"Unknown symmetrization {symmetrize!r}.")

    def decorator(body: Callable) -> Callable:
        body_signature = inspect.signature(body)
        body_params = list(body_signature.parameters.values())
        combine = SYMMETRIZE[symmetrize] if symmetrize else None
        similarity = kind == "simil"

        @functools.wraps(body)
        def wrapper(seq_x, seq_y, **kwargs):
            use_normal = kwargs.pop("normal", False) if normal else False
            use_directional = (
                kwargs.pop("directional", False) if directional_option else False
            )

            # Validate the parameters first, applying their defaults
            bound_args = body_signature.bind(seq_x, seq_y, **kwargs)
            bound_args.apply_defaults()
            options = dict(list(bound_args.arguments.items())[2:])
            if check is not None:
                check(**options)

            # Empty rule
            if empty == "max":
                value = empty_dissim(seq_x, seq_y)
                if value is not None:
                    return 1.0 - value if similarity else value

            if bound == "custom":
                options["normal"] = use_normal

            # Raw value, symmetrized if needed
            raw = body(seq_x, seq_y, **options)
            if combine is not None and not use_directional:
                other = body(seq_y, seq_x, **options)
                if isinstance(raw, Scored):
                    raw = Scored(combine(raw.value, other.value), raw.bound)
                else:
                    raw = combine(raw, other)

            # Normalization and float result
            if bound == "scored":
                if use_normal:
                    return raw.value / raw.bound if raw.bound else 0.0
                return float(raw.value)
            if use_normal and bound in ("max_len", "sum_len"):
                size = (
                    max(len(seq_x), len(seq_y))
                    if bound == "max_len"
                    else len(seq_x) + len(seq_y)
                )
                return raw / size if size else 0.0
            if use_normal and bound == "clip":
                return min(max(raw, 0.0), 1.0)

            return float(raw)

        # The public signature: the parameters of the body, plus
        # `directional` and `normal`
        params = [
            param
            for param in body_params
            if not (bound == "custom" and param.name == "normal")
        ]
        keyword = inspect.Parameter.KEYWORD_ONLY
        if directional_option:
            params.append(
                inspect.Parameter(
                    "directional", keyword, default=False, annotation=bool
                )
            )
        if normal:
            params.append(
                inspect.Parameter("normal", keyword, default=False, annotation=bool)
            )
        wrapper.__signature__ = body_signature.replace(parameters=params)

        if normal and body.__doc__:
            wrapper.__doc__ = _document_normal(
                body.__doc__, _describe_normal(bound, normal_doc)
            )

        module = body.__module__.rsplit(".", 1)[-1]
        info = MeasureInfo(
            func=wrapper,
            name=f"{module}.{body.__name__}",
            key=key,
            kind=kind,
            identity=identity if kind in ("dist", "dissim") else None,
            identity_example=identity_example,
            identical_zero=identical_zero,
            symmetric=symmetric,
            triangle=triangle if kind in ("dist", "dissim") else None,
            triangle_example=triangle_example,
            condition=condition,
            bound=bound,
            raw_range=raw_range or _default_range(bound),
            empty=empty,
            normal_effect=normal and bound not in ("unit",),
            accepts_normal=normal,
        )
        wrapper.info = info
        _REGISTRY.append(info)

        return wrapper

    return decorator


def _default_range(bound: str) -> str:
    """
    Returns a description of the range of raw values for a bound.
    """

    return {
        "max_len": "0 to max length",
        "sum_len": "0 to sum of lengths",
        "unit": "0 to 1",
        "clip": "about 0 to 1",
    }.get(bound, "0 upwards")


def measures(
    kind: Optional[str] = None, *, in_distance: Optional[bool] = None
) -> List[MeasureInfo]:
    """
    Returns the declarations of the measures of the library.

    Each declaration (a `MeasureInfo`) holds the function, its name in
    `distance()` (if any), its kind, and its claims: identity of
    indiscernibles, symmetry, and triangle inequality, with counterexamples
    for the claims that fail.

    Example
    ********

    .. code-block:: python

        >>> [info.key for info in seqsim.measures("dist")][:3]
        ['levenshtein', 'levenshtein_gld', 'levenshtein_ned']

    :param kind: Only return measures of this kind (see `KINDS`).
    :param in_distance: If given, only return measures available (`True`)
        or not available (`False`) through `distance()`.
    :return: The declarations, in order of registration.
    """

    # Make sure all modules with measures are imported
    from . import alignment, compression, edit, order, sequence, token  # noqa: F401

    return [
        info
        for info in _REGISTRY
        if (kind is None or info.kind == kind)
        and (in_distance is None or (info.key is not None) == in_distance)
    ]


def methods() -> Dict[str, Callable]:
    """
    Returns the mapping of method names to functions, sorted by name.
    """

    return {
        info.key: info.func
        for info in sorted(measures(in_distance=True), key=lambda info: info.key)
    }


def _cell(info: MeasureInfo, field: str) -> str:
    """
    Returns the text of a cell of the table of measures.
    """

    claim = info.identity if field == "identity" else info.triangle
    if claim == "no":
        example = (
            info.identity_example if field == "identity" else info.triangle_example
        )
        return "no (e.g. " + ", ".join(f"`{seq!r}`" for seq in example) + ")"
    if claim == "conditional":
        return info.condition

    return {"unproven": "not proven"}.get(claim, claim)


def markdown_table(link: Optional[Callable[[MeasureInfo], str]] = None) -> str:
    """
    Returns the Markdown table of the measures available through `distance()`.

    :param link: An optional function returning the link target of the name
        of each function.
    :return: The table, one line per measure, with a header.
    """

    lines = [
        "| `distance()` key | Function | Identity of indiscernibles | "
        "Triangle inequality | Raw range |",
        "|---|---|---|---|---|",
    ]
    ordered = sorted(
        measures(in_distance=True),
        key=lambda info: (KINDS.index(info.kind), info.name),
    )
    for info in ordered:
        name = f"`{info.name}`"
        if link is not None:
            name = f"[{name}]({link(info)})"
        lines.append(
            f"| `{info.key}` | {name} | {_cell(info, 'identity')} | "
            f"{_cell(info, 'triangle')} | {info.raw_range} |"
        )

    return "\n".join(lines) + "\n"


__all__ = ["measure", "measures", "methods", "markdown_table", "MeasureInfo", "Scored"]
