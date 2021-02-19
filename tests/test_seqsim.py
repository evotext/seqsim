#!/usr/bin/env python

"""
test_seqsim
===========

Tests for the `seqsim` module.
"""

import unittest
from src import seqsim
from src.seqsim import Simhash

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

    # simhash testing
    def test_simhash_distance(self):
        sh = Simhash('How are you? I AM fine. Thanks. And you?')
        sh2 = Simhash('How old are you ? :-) i am fine. Thanks. And you?')
        self.assertTrue(sh.distance(sh2) > 0)

        sh3 = Simhash(sh2)
        self.assertEqual(sh2.distance(sh3), 0)
        self.assertNotEqual(Simhash('1').distance(Simhash('2')), 0)

    def test_simhash_chinese(self):
        self.maxDiff = None

        sh1 = Simhash(u'你好　世界！　　呼噜。')
        sh2 = Simhash(u'你好，世界　呼噜')

        sh4 = Simhash(u'How are you? I Am fine. ablar ablar xyz blar blar blar  blar blar blar blar Thanks.')
        sh5 = Simhash(u'How are you i am fine.ablar ablar xyz blar blar blar blar blar blar blar than')
        sh6 = Simhash(u'How are you i am fine.ablar ablar xyz blar blar blar blar blar blar ablar thank')

        self.assertEqual(sh1.distance(sh2), 0)

        self.assertTrue(sh4.distance(sh6) < 3)
        self.assertTrue(sh5.distance(sh6) < 3)

COMMENTED = """

    def test_hayneedle(self):
        # read data from the Unix dictionary
        with open('/usr/share/dict/words') as handler:
            words = [line.strip() for line in handler]

        # pop some random needles
        random.shuffle(words)
        CUT = 2
        needles, hay = words[:CUT], words[CUT:]

        for method in ['mmcwpa', 'simhash']:
            dists = seqsim.stringscomp(needles, hay, method=method)
            for needle in needles:
                print(method, needle, dists[needle][:3])
"""


if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())

