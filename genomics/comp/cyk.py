import numpy as np
import nltk

"""
Cocke-Younger-Kasami (CYK Cocke younger kasami) algorithm to calculate
the most optimal parse tree for a sequence. Create a Well-Formed 
Substring Table (WFST wfst).
"""

tokens = ["the", "kids", "opened", "the", "box", "on", "the", "floor"]

grammar = nltk.parse_cfg("""
    S -> NP VP
    PP -> P NP
    NP -> Det N | NP PP
    VP -> V NP | VP PP
    Det -> "the"
    N -> "kids" | "box" | "floor"
    V -> "opened"
    P -> "on"
    """)

def init_wfst(tokens, grammar):
    numtokens = len(tokens)
    # fill with dots
    wfst = [["." for i in range(numtokens+1)]
                 for j in range(numtokens+1)]
    # fill diagonal
    for i in range(numtokens):
        productions = grammar.productions(rhs=tokens[i])
        wfst[i][i+1] = productions[0].lhs()
    return wfst 
