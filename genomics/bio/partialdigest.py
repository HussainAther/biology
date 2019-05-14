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

def spd(l):
    """
    For some length list of dna l, find the actual places of cut.
    """
