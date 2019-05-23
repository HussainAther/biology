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
    """
    Initialize the WFST table.
    """
    numtokens = len(tokens)
    # fill with dots
    wfst = [["." for i in range(numtokens+1)]
                 for j in range(numtokens+1)]
    # fill diagonal
    for i in range(numtokens):
        productions = grammar.productions(rhs=tokens[i])
        wfst[i][i+1] = productions[0].lhs()
    return wfst 

def complete_wfst(wfst, tokens, trace=False):
    """
    Fill in the rest of the table.
    """
    index = {}
    for prod in grammar.productions(): # make reverse lookup
        index[prod.rhs()] = prod.lhs()
    numtokens = len(tokens)
    for span in range(2, numtokens+1):
        for start in range(numtokens+1-span): # go down diagonal
            end = start + span
            for mid in range(start+1, end):
                nt1, nt2 = wfst[start][mid], wfst[mid][end]
                if (nt1,nt2) in index:
                    if trace:
                        print("[%s] %3s [%s] %3s [%s] ==> [%s] %3s [%s]" % \
                               (start, nt1, mid, nt2, end, start, index[(nt1,nt2)], end))
                        wfst[start][end] = index[(nt1,nt2)]
    return wfst 
