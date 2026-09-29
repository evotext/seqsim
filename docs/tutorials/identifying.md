# Identifying texts across manuscripts

Before the contents of manuscripts can be compared, the texts they contain
must be identified: the same text appears under different titles, rubrics,
and incipits, with variable spelling, abbreviations, and word order. This
tutorial uses `seqsim` to propose matches between the texts of two
manuscripts, a task that usually precedes the analyses of the
[contents tutorial](contents.md).

## The data

The rubrics of two invented manuscripts of saints' lives:

```python
>>> import seqsim
>>> ms_a = [
...     "Incipit vita sancti Antonii abbatis",
...     "Vita beati Pauli primi eremite",
...     "Passio sanctae Agnetis virginis et martyris",
...     "Vita sancti Martini episcopi",
... ]
>>> ms_b = [
...     "Vita s. Martini Turonensis episcopi",
...     "Incipit passio beate Agnetis uirginis",
...     "Vita Pauli heremite",
...     "Vita et conversatio beati Antonii",
... ]
```

## Comparing single rubrics

Rubrics can be compared as strings of characters. Normalizing the case and
the spelling (here, *u*/*v* and *e*/*ae*, and the removal of *h*) first
avoids counting trivial differences:

```python
>>> def normalize(text):
...     text = text.lower().replace("v", "u").replace("ae", "e").replace("h", "")
...     return text
>>> seqsim.edit.levenshtein_dist(normalize(ms_a[1]), normalize(ms_b[2]), normal=True)
0.4
>>> seqsim.edit.levenshtein_dist(normalize(ms_a[1]), normalize(ms_b[0]), normal=True)
0.6857142857142857
```

Rubrics differ, however, in words added or omitted (*sancti*, *beati*,
*incipit*, epithets), and in word order, which character-level measures
handle poorly.

## Comparing rubrics word by word

The Monge-Elkan similarity treats each rubric as a sequence of words: for each
word of one rubric, it finds the most similar word of the other, and averages
these values. It tolerates different word orders and, with a character-level
measure for the words, spelling differences.

```python
>>> def words(text):
...     return normalize(text).split()
>>> seqsim.alignment.monge_elkan_simil(words(ms_a[0]), words(ms_b[3]))
0.6284632034632035
>>> seqsim.alignment.monge_elkan_simil(words(ms_a[0]), words(ms_b[0]))
0.45595238095238094
```

Common words, such as *vita*, *incipit*, *sancti*, or *beati*, carry little
information and can make unrelated rubrics look similar. Removing them before
the comparison focuses on the distinctive words:

```python
>>> stopwords = {"incipit", "uita", "passio", "sancti", "sancte", "beati",
...              "beate", "s.", "et", "conuersatio"}
>>> def keywords(text):
...     return [word for word in words(text) if word not in stopwords]
>>> keywords(ms_a[2])
['agnetis', 'uirginis', 'martyris']
>>> seqsim.alignment.monge_elkan_simil(keywords(ms_a[2]), keywords(ms_b[1]))
0.8958333333333333
```

## Proposing matches

With a similarity for pairs of rubrics, each text of one manuscript can be
matched to its most similar text in the other:

```python
>>> def similarity(rubric_a, rubric_b):
...     return seqsim.alignment.monge_elkan_simil(keywords(rubric_a), keywords(rubric_b))
>>> for rubric in ms_a:
...     best = max(ms_b, key=lambda other: similarity(rubric, other))
...     print(f"{rubric[:30]:<30} -> {best[:30]:<30} {similarity(rubric, best):.2f}")
Incipit vita sancti Antonii ab -> Vita et conversatio beati Anto 0.82
Vita beati Pauli primi eremite -> Vita Pauli heremite            0.90
Passio sanctae Agnetis virgini -> Incipit passio beate Agnetis u 0.90
Vita sancti Martini episcopi   -> Vita s. Martini Turonensis epi 0.88
```

All four texts are correctly identified. In real collections, such
automatic matches should be treated as proposals to be checked, for example
by reviewing all matches below a threshold, or where the best and second best
candidates are close.

## Other measures

- The Jaro-Winkler dissimilarity (`edit.jaro_winkler_dissim`) favours strings
  sharing their beginning, and is well suited to comparing single names, such
  as the name of a saint or an author.
- The global alignment (`alignment.nw_dissim`) allows custom costs, for
  example to make the substitution of letters commonly confused by scribes
  or readers (*c*/*t*, *u*/*n*, *f*/*s*) cheaper.
- Incipits can be compared as sequences of words with the edit distances,
  as in the [witnesses tutorial](witnesses.md), since their word order is
  usually stable.
