#!/usr/bin/env python

"""
test_seqsim
===========

Tests for the `seqsim` module.
"""

import unittest
import seqsim

class TestSeqsim(unittest.TestCase):

#    tests = [
#        ['abc', 'def'],
#        ['abcdef', 'abcdef'],
#        ['Austria', 'Australia'],
#        ['Python', 'python'],
#        ['a123b', 'ab123'],
#        ['129 Industry Park', '129 Indisttry Park'],
#        ['abc de', 'abc k de'],
#        ['de abc', 'de abc'],
#        ['Fu Hui', 'Mr Fu Hui'],
#        ['Fu Hui', 'Fu Mr Hui'],
#        ['abcdefgh ijklmnpo', 'abcdefh ijklmnwo'],
#        ['Gao Hua Ming', 'Gao Ming Hua'],
#        ['zeng zeng', 'zeng hong'],
#    ]

    def test_000_seqsim(self):
        self.assertTrue(seqsim.stringcomp('abc', 'def') == 0.0)

    def test_001_seqsim(self):
        self.assertTrue(seqsim.stringcomp('abc', 'abc') == 1.0)

    def test_002_seqsim(self):
        self.assertTrue(seqsim.stringcomp('Python', 'python') > 0.5)


if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())

