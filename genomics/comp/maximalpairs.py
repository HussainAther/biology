import numpy as np

"""
Find all-right maximal pairs in string s with bounded gap.
"""

def report(a, b):
    """
    Return the  
    """
    result = []

def pairs():
    """
    1. Build the binary suffix tree T and create at each leaf an AVL tree of size one
       that stores the index at the leaf. An AVL tree is a blaanced search tree that stores
       an ordered set of elements. 
    2. When the AVL trees T1 and T2 (in which |T1| <= |T2|) at the two children w1 and w2 of
       a node v with path-label alpha are available we do:
       (a) Let {p1, p2,...} be the elements in T1 in sorted order. For each element p in T1, we find
           qr(p) = min{x of set T2 | x >= p + |alpha| + g1(|alpha|)}
           ql(p) = min{x of set T2 | x >= p - |alpha| - g2(|alpha|)}
           after searching in T2 with the two sorted lists {pi + |alpha| + g1(|alpha|) i = 
           {1, 2,..s} and {pi - |alpha| - g2(|alpha|) | i = 1, 2,...s}
       (b) For each element p in T2 we call report(qr(p), p + |alpha| + g2(|alpha|)) and
           report (as defined above).
       (c) Build the leaf-list tree T at node v by merging T1 and T2 applying the following lemma:
          "Two AVL trees of size at most n and m, where n <= m, can be merged together
           in time O(nlog(m/n))"
    """
