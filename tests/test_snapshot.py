"""
test_snapshot
=============

Checks that every public measure returns exactly the recorded outputs (see
`snapshot_cases.py`), so that refactoring cannot change any result.
"""

# Import Python standard libraries
import json

# Import the cases
from snapshot_cases import SNAPSHOT, compute


def test_snapshot():
    expected = json.loads(SNAPSHOT.read_text())
    actual = compute()

    assert sorted(actual) == sorted(expected)
    differences = {
        key: (expected[key], actual[key])
        for key in expected
        if actual[key] != expected[key]
    }
    assert (
        not differences
    ), f"{len(differences)} changed results, e.g. {list(differences.items())[:5]}"
