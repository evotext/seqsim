# encoding: utf-8

__author__ = "Tiago Tresoldi"
__email__ = "tresoldi@shh.mpg.de"
__version__ = "0.2"

from .mmcwpa import mmcwpa
from .simhash import Simhash, SimhashIndex

def stringcomp(str_x, str_y, method='mmcwpa'):
    if method == 'mmcwpa':
        ret = mmcwpa(str_x, str_y)
    elif method == 'simhash':
        print("Simhash implementaiton missing")
    else:
        print("Error, invalid method")

    return ret

def stringscomp(needles, hay, method='simhash'):
    if method == 'mmcwpa':
        # brute force compare the hay with all needles and return
        dists = {}
        for needle in needles:
            dists[needle] = {seq:stringcomp(needle, seq) for seq in hay}

        # get the sorted ones
        sorted_dists = {
            needle : sorted(dists[needle].items(), key=lambda k:k[1], reverse=True)
            for needle in dists}

    elif method == 'simhash':
        # compute all hashes
        hashes = {word : Simhash(word) for word in hay}

        dists = {}
        for needle in needles:
            needle_hash = Simhash(needle)
            dists[needle] = {
                seq : needle_hash.distance(hashes[seq]) for seq in hay
            }

        sorted_dists = {
            needle : sorted(dists[needle].items(), key=lambda k:k[1])
            for needle in dists
        }
    else:
        print("error method not implemented")

    return sorted_dists
