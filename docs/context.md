# seqsim: domain glossary

The vocabulary used in the code, the documentation, and design discussions.
Module and file names follow these terms.

## Comparing two sequences

**Sequence**
: An ordered collection of elements (a string, list, or tuple). The library
  never looks inside an element.

**Element / item**
: A hashable Python object in a sequence. Two elements are the same if they
  are equal (`==`). "Item" is preferred when the elements are units of
  content, such as the texts of a manuscript.

**Occurrence**
: A repeated element is distinguished by its occurrence: the first `"a"`, the
  second `"a"`, and so on. Measures of order match occurrences in order
  (`repeats="occurrence"`); some analyses keep only the first occurrence of
  each element (`repeats="first"`).

**Adjacency**
: A pair of consecutive elements `(x, y)`. With **boundaries**, a start and
  an end marker are added, so that a sequence of `n` elements has `n + 1`
  adjacencies.

**n-gram (shingle, q-gram)**
: A contiguous sub-sequence of `n` elements, optionally padded with a
  boundary marker at both ends.

**Measure**
: A function comparing two sequences, `f(x, y, *, normal=False, ...)`,
  returning a float. Every measure is registered with its metadata (see
  `seqsim.measures()`).

**Kind**
: What a measure returns, stated by the suffix of its name:
  - `_dist`, a **distance**: a true metric (non-negativity, identity of
    indiscernibles, symmetry, triangle inequality);
  - `_dissim`, a **dissimilarity**: 0 for identical sequences, higher for
    more different ones, metric properties not all guaranteed;
  - `_simil`, a **similarity**: higher for more similar sequences.
  Functions outside this scheme are **directional** (such as
  `token.containment`) or **estimators** (such as `order.iebp_estimate`).

**Claim**
: A mathematical property a measure is declared to have or to lack
  (identity of indiscernibles, triangle inequality, whether `normal` has an
  effect). Claims that hold are tested; claims that fail carry a
  **counterexample**, which is tested too.

**Raw value / normalized value**
: A measure returns its raw value by default, and with `normal=True` a value
  in range [0..1], obtained by dividing by its **bound**, an upper limit of
  the raw value for the pair of sequences (such as the length of the longest
  sequence).

**Symmetrization**
: Measures whose algorithm depends on the order of the arguments are
  computed in both orders, keeping the result indicating the greatest
  similarity, so that the measure is symmetric.

**Empty rule**
: How a measure treats empty sequences: two empty sequences are identical,
  and an empty sequence compared with a non-empty one has the largest value
  (for dissimilarities in range [0..1], `1.0`), unless the natural value of
  the algorithm (such as the cost of inserting every element) applies.

## Traditions

**Witness**
: A sequence of items representing one manuscript, print, or version, such
  as the list of texts in a codex.

**Tradition**
: A collection of witnesses of the same work or collection, mapped from
  witness names to witnesses.

**Part**
: A section of a codex catalogued separately; complementary parts of one
  codex can be merged into a single witness.

**Frame**
: The reference order of all the items of a tradition, such as a consensus
  order.

**Label**
: A unit of the frame (a chapter, a book) assigned to each item.

**Coverage**
: For each item of the frame and one witness, whether the witness has it
  (`1`), lacks it where its absence is evidence (`0`), or is missing data
  (`None`), for example because of a lacuna.

**Character**
: A column of a character matrix for phylogenetic software: a **content**
  character (whether a witness has an item) or an **adjacency** character
  (whether an item is immediately followed by another).
