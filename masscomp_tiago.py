import sys
from src import seqsim

import random

def test_hayneedle(method, num_needles):
    # read data from the Unix dictionary
    with open('/usr/share/dict/words') as handler:
        words = [line.strip() for line in handler]

    # pop some random needles
    random.shuffle(words)
    needles, hay = words[:num_needles], words[num_needles:]

    # run
    dists = seqsim.stringscomp(needles, hay, method=method)
    for needle in needles:
        print(method, needle, dists[needle][:3])

if __name__ == "__main__":
    method = sys.argv[1]
    needles = int(sys.argv[2])

    test_hayneedle(method, needles)
