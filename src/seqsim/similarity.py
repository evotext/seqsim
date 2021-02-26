import difflib


def birnbaum(M, N):
    sm = difflib.SequenceMatcher(None, M, N)
    blocks = sm.get_matching_blocks()
    sizes = [1 if match.size == 0 else match.size for match in blocks]
    formula = [(v * (v + 1)) // 2 for v in sizes]

    return sum(formula)
