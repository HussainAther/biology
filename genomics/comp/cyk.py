import numpy as np
import nltk

"""
Cocke-Younger-Kasami (CYK Cocke younger kasami) algorithm to calculate
the most optimal parse tree for a sequence.
"""

tokens = ["the", "kids", "opened", "the", "box", "on", "the", "floor"]

grammar = nltk.parse_cfg("""
    S -> NP VP
    PP -> P NP
    NP -> Det N | NP PP
    VP -> V NP | VP PP
    Det -> 'the'
    N -> 'kids' | 'box' | 'floor'
    V -> 'opened'
    P -> 'on'
    """)
