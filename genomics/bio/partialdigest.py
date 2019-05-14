"""
The Simplified Partial Digest Problem (simplified partial digest) uses primary fragments and base fragments to locate 
restriction sites. Base fragments have two endpoints that were consecutive sites on the
target DNA strand and can be obtained by exposing the strand to the enzyme until the digestion process is
complete. We consider there are n restriction sites where the enzyme cuts along a DNA strand of length D.
Simplified Partial Digest Problem (SPDP) Statement: Given X0 = 0, Xn+1 = D, and a set of base
fragments {Xi - Xi-1}1 ≤ i ≤ n+1 and primary fragments {(Xn+1 - Xi) ∪ Xi}1 ≤ i ≤ n, reconstruct the original series
X1,...,Xn, where Xi corresponds to the distance between the leftmost end of the target DNA strand and the ith
furthest restriction site along the strand. 
"""

def delete(elements, a): 
    """
    Remove elements from list a.
    """
    for el in elements:
        a.remove(el)
    return a

def delta(y, X):
    """
    For X list and integer y, subtract, find the absolute values, and sort.
    """ 
    n = len(X)
    for i in range(n):
        X[i] -= y
        X[i] = abs(X[i])
    return sorted(X)

def contains(small, big):
    """
    If small has big.
    """
    for i in range(len(big) - len(small)+1):
        for j in range(len(small)):
            if big[i+j] != small[j]:
                break
        else:
            return True
    return False

def partialDigest(L):
    """
    Partial digest algorithm.
    """
    global width
    width = (max(L)) 
    delete([width], L) # Needs to be in list to feed to 'delete' function
    X = [0, width]
    X = place(L,X)
    return X

def place(L, X):
    """
    Place L in X.
    """
    if len(L) == 0: # Baseline condition
        return X
    y = max(L)

    if contains(delta(y,X),L): # If former is the subset of L
        delete(delta(y,X), L)  # Remove lengths from L
        X += list(y) # assert that this y is one of the fixed points, X
        X = sorted(X) # To maintain order
        print(X)

        place(L,X) # Recursive call of the function to redo the upper part
                       # If none of the if statements match the condition, continue

        X.remove(y) # If the code reaches down here, it means the assumption that
                        # y is one of the points is wrong. Thus undo
        L += delta(y,X) # undo L
        L = sorted(L)   # To maintain order

    # Do the same thing except this time it's (width-y)
    elif contains(delta(width-y,X),L): 
        delete(delta(width-y,X), L)
        X += list(width - y) 
        X = sorted(X) 

        place(L,X) 

        X.remove(width-y) 

        L += delta(y,X) 
        L = sorted(L)
    return X

