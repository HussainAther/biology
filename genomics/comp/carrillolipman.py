"""
Carrillo–Lipman multiple alignment algorithm (carrillo lipman).

Substitution costs.
S(A,A) = S(C,C) = S(G,G) = S(T,T) = 0, 
S(C,T) = S(T,C) = S(A,G) = S(G,A) = 1,
S(A,T) = S(T,A) = S(A,C) = S(C,A)
       = S(G,C) = S(C,G) = S(G,T) = S(T,G) = 2

with gap-open cost d = 3 and gap-extension cost e = 2.
"""
