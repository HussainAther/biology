import unittest

from dnaseqlib import *

"""
DNA Matching utility functions.
"""

class RollingHash:
    """
    Utility class
    """
    def __init__(self, s):
        """
        Initialize self with a given sequence s.
        """
        self.HASH_BASE = 7
        self.seqlen = len(s)
        n = self.seqlen - 1
        h = 0
        for c in s:
            h += ord(c) * (self.HASH_BASE ** n)
            n -= 1
        self.curhash = h

    def hash(self):
        """
        Return the current hash value.
        """
        return self.curhash

    def slide(self, previtm, nextitm):
        """
        Update the hash by removing the previous item previtm and adding the next item nextitm.
        Return the updateed hash value.
        """
        self.curhash = (self.curhash * self.HASH_BASE) + ord(nextitm)
        self.curhash -= ord(previtm) * (self.HASH_BASE ** self.seqlen)
        return self.curhash

class Multidict:
    """
    Create a multidict from several seuqences.
    """
    def __init__(self, pairs=[]):
        """
        Initialize the self with a list of pairs.
        """
        self.table = dict()
        for pair in pairs:
            self.put(pair[0], pair[1])

    def put(self, k, v):
        """
        Put key k and value v in the table or append it to an existing entry.
        """
        if k in self.table:
            self.table[k].append(v)
        else:
            self.table[k] = [v]

    def get(self, k):
        """
        Retrieve a key k from the table.
        """
        try:
            return self.table[k]
        except KeyError:
            return []

"""
Searches for commonalities between sequences a and b by comparing subsequences
of length k.  The sequences a and b should be iterators that return nucleotides.
The table is built by computing one hash every m nucleotides (for m >= k).
"""

def getExactSubmatches(a, b, k, m):
    """
    Get the submatches from the multidict adding the key k with a value m
    for the multidict a. Create hashes of the sequence b for a key k
    such that we can compare. 
    """
    print("Building table from sequence A...")
    seqtable = Multidict(intervalSubsequenceHashes(a, k, m))
    print("...done building table.")
    for hashval, (bpos, bsubseq) in subsequenceHashes(b, k):
        for apos, asubseq in seqtable.get(hashval):
            if asubseq != bsubseq:
                continue
            yield (apos, bpos)
    return

def subsequenceHashes(seq, k):
    """
    Generate subsequence hashes for comparing subsequences for 
    a key k and sequence seq.
    """
    try:
        assert k > 0
        subseq = ""
        for i in range(0, k):
            subseq += seq.next()
        rh = RollingHash(subseq)
        pos = 0
        while True:
            yield (rh.hash(), (pos, subseq))
            previtm = subseq[0]
            subseq = subseq[1:] + seq.next()
            rh.slide(previtm, subseq[-1:])
            pos += 1
    except StopIteration:
        return

def intervalSubsequenceHashes(seq, k, m):
    """
    Use the interval-method to manage subsequence hashes.
    """
    assert m >= k
    try:
        pos = 0
        while True:
            subseq = ""
            for i in range(0, k):
                subseq += seq.next()
            rh = RollingHash(subseq)
#            print (rh.hash(), (pos, subseq))
            yield (rh.hash(), (pos, subseq))
            for i in range(0, m-k):
                seq.next()
            pos += m
    except StopIteration:
        return

# Run this to return the exact matches between the sequences.
#compareSequences(getExactSubmatches, "mouse-paternal_8_100.png", (500,500), "fmouse.fa", "fpaternal.fa", 8, 100)

class TestRollingHash(unittest.TestCase):
    """
    Create a rolling hash to test the exact matches between the sequences.
    """
    def test_rolling(self):
        rh1 = RollingHash("abcde")
        rh2 = RollingHash("bcdef")
        rh3 = RollingHash("cdefZ")
        rh1.slide("a","f")
        self.assertTrue(rh1.hash() == rh2.hash())
        rh1.slide("b","Z")
        self.assertTrue(rh1.hash() == rh3.hash())

class TestMultidict(unittest.TestCase):
    def test_multi(self):
        foo = Multidict()
        foo.put(1, "a")
        foo.put(2, "b")
        foo.put(1, "c")
        self.assertTrue(foo.get(1) == ["a","c"])
        self.assertTrue(foo.get(2) == ["b"])
        self.assertTrue(foo.get(3) == [])

if __name__ == "__main__":
    unittest.main()
