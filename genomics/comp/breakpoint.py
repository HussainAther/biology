"""
We may use breakpoints when dealing with permutations to figure out
what the identity permutation (ancestor) from which they originated.
We assume the ancestor doesn't have any breakpoints and note the 
breakpoints that appear as reversals occur in the sequence. 
"""

def isSorted(l):
    """
    Is the list l sorted?
    """
    return sorted(l) == list(range(min(l), max(l)+1))

def hasDecrease(l):
    """
    Is there a decreasing strip of numbers in l?
    """
    for i in range(len(l) - 1):
        if l[i+1] < l[i]:
            return True
    return False

def reversal(l):
    """
    Apply the reversal step throughout l to remove most breakpoints.
    """
    bpi = "" # breakpoint initial index
    bpf = "" # breakpoint final index
    bpif = [] # list of tuples of initial and final indices of breakpoints for reversals 
    for i, j in enumerate(l):
        if l[i+1] != l[i] + 1 and bpi == "":
            bpi = i
        elif l[i+1] != l[i] + 1 and bpf != "":
            bpf = i
        elif bpi != "" and bpf != "":
            bpif.append((bpi, bpf))
    for i, f in bpif:
        l[i:f] = l[i:f][::-1] # reverse the strip 
             
def reverseLongestIncrease(l):
    """
    Reverse the longest increasing strip in l.
    """
    strip = [] # indices of the strip
    tempstrip = [] # indices of the temporary strip used for comparison
    for i in range(len(l) - 1):
        if l[i+1] > l[i]:
            tempstrip.append(i)
        elif tempstrip != []:
            if len(tempstrip) > len(strip):
                strip = tempstrip
    l[strip[0]:strip[-1]] = l[strip[0]:strip[-1]][::-1] # reverse the strip 

def brs(a):
    """
    For an array a of sequences of numbers, remove as many breakpoints as possible
    in every step until we get the numbers in order. 
    """
    for i in a:
        while not isSorted(i):
            if hasDecrase(i):
                reversal(i)
            else:
                reverseLongestIncrease(i)
     
