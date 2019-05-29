"""
Sequence by hybridization (SBH) sequencing approach that uses a miniature DNA array
of a DNA chip that has thousands of short DNA fragments called probes. Given a short 
probe (an 8- to 30-nucleotide single-stranded synthetic DNA fragment) and a single-stranded 
target DNA fragment, the target will hybridize with the probe if the probe is a substring 
of the target’s Watson-Crick complement. When the probe and the target are mixed together, 
they form a weak chemical bond and stick together. For example, a probe ACCGTGGA
will hybridize to CCCTGGCACCTA since it's complementary to the substring TGGCACCT of the target.

For a string s of length n, the l-mer composition, or spectrum, of s, is the multiset of n − l + 1 l-mers 
in s and is written Spectrum(s, l). If l = 3 and s = TATGGTGC, then Spectrum(s, l) = {TAT, ATG, TGG, GGT, GTG, TGC}. 
We can now formulate the problem of sequencing a target DNA fragment from its DNA array data.
"""

def superstring(a, b): 
    """
    Return the shortest common supersequence (superstring) between a and b.
    """
    m = len(a) 
    n = len(b) 
    l = lcs(a, b, m, n) 
    # Result is sum of input string 
    # lengths - length of lcs 
    return (m + n - l)

def lcs(a, b, m, n):
    """
    Return lcs length for a and b of lengths m and n respectively.
    """ 
    L = [[0] * (n + 2) for i in
                    range(m + 2)] 
    # Following steps build L[m + 1][n + 1] 
    # in bottom up fashion. Note that L[i][j] 
    # contains length of LCS of a[0..i - 1] 
    # and b[0..j - 1] 
    for i in range(m + 1): 
        for j in range(n + 1): 
            if (i == 0 or j == 0) : L[i][j] = 0
            elif (a[i - 1] == b[j - 1]) : 
                L[i][j] = L[i - 1][j - 1] + 1
            else : L[i][j] = max(L[i - 1][j], 
                                 L[i][j - 1]) 
    # L[m][n] contains length of 
    # LCS for a[0..n - 1] and b[0..m - 1] 
    return L[m][n]  

def sbh(x, l):
    """
    For a list x of all l-mers of a string s, return string s such that Spectrum(s, l) = x. 
    """
    s = ""
    t = iter(x)
    for i in t:
        print(i, next(t))
