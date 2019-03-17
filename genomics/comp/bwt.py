# -*- coding: utf-8 -*-
#!/usr/bin/python
import os, re, sys

from itertools import islice

"""
Burrows-Wheeler alignment transform.
"""

class BurrowsWheeler():
    EOS = "\0" #please use this marker, it's the "null" character, and has been carefully chosen to sort before all other characters
    # EOS = "#" # a visible end marker
    
    def transform(self, s):
        """ Simplest Burrows-Wheeler transform implementation, O(n^2) respective
            to the length of the text. """
        assert self.EOS not in s, "Input string cannot contain null character (%s)" % self.EOS
        
        # add end of text marker
        #please note that although this code uses the null marker to mark the end of the string,
        #i've replaced it with "$" just for the output to help you debug your code -BZ
        s += self.EOS
        
        # s is the input string, with the end-of-string character added for you
        # generate a list of rotated input strings
        rotations=[]
        for i in range(len(s)):
            rotations.append(s[i:]+s[:i])
            # sort the table of rotations
            sortedRotations=sorted(rotations)
        
        # extract the last characters of each row
        # convert the characters to a string, and store them in the variable bwtTransformedString
        bwtTransformedString=''
        for rotation in sortedRotations:
            bwtTransformedString+=rotation[-1]
        return bwtTransformedString
    
    def inverse(self, s):
        """ Simplest Inverse Burrow-Wheeler transform implementation. """
        # make empty table for the suffix array
        table = [""] * len(s)
        
        # use lf-mapping to reverse the tranformation
        for i in range(len(s)):
            # add one letter for each partial string in the suffix array
            prepended = [s[i] + table[i] for i in range(len(s))]
            
            # convert it to sorted suffix array
            table = sorted(prepended)
            
        # Find the correct row (ending in "\0")
        for row in table:
            if row.endswith(self.EOS):
                s = row
                break
        
        # Get rid of trailing null character
        s = s.rstrip(self.EOS)
        
        return s
# ---------------------------------------------------------------------------- #
# Different Inverse implementations
# ---------------------------------------------------------------------------- #


def calc_first_occ(s):
    #returns C[c], a dictionary for each letter c, of the number of lexically smaller letters in the string -BZ
    """ calculate the first occurance of a letter in sorted string s """
    # s - is the bwt transformed string -or the original, it doesn't matter
    A = {} # letter count
    #get the number of times each letter occurs
    for i, c in enumerate(s): #enumerate returns a list of enumerated tuples eg: enumerate("ABC")-->[(0,'A'),(1,'B'), (2,'C')]
        if A.get(c):
            A[c] += 1
        else:
            A[c] = 1
    
    # sort the letters
    letters = sorted(A.keys())
    
    # first index of letter
    occ = {}
    
    #iterate through letters in order, incrementing index as we go, to get index of first occurrence of each letter.
    idx = 0
    for c in letters:
        occ[c] = idx
        idx += A[c]
    del idx, A
    
    return occ

def calc_checkpoints(s, step):
    #returns the list of checkpoints, seperated by distance "step", each of which is a dictionary of the
    #number of occurences of each character up to that point (noninclusive) -BZ
    """ count the number of letters for each step and
        return list of the counts"""
    A = {} # letter count
    C = [] # checkpoints
    for i, c in enumerate(s):
        if i % step == 0:
            C.append(A.copy())
        if A.get(c):
            A[c] += 1
        else:
            A[c] = 1
    return C

def count_letter_with_checkpoints(C, step, s, idx, letter):
    """ Count the number of a letter upto idx in s using checkpoints.
    
    Arguments:
    C      -- is the list of checkpoints
    step   -- is the step of the checkpoints
    s      -- the transformed string
    idx    -- count upto this position
    letter -- count for this letter
    """
    
    # find the nearest checkpoint for idx
    check = (idx + (step / 2)) / step
    if check >= len(C):
        check = len(C) - 1
    pos = check * step
    
    # count of the letter s[idx] upto pos (not included)
    count = C[check].get(letter)
    if count == None:
        count = 0
    
    # range between pos and idx
    if pos < idx:
        r = xrange(pos, idx)
    else:
        r = xrange(idx, pos)
    
    # count of letters between pos, idx
    k = 0        
    for i in r:
        if letter == s[i]:
            k += 1
    
    # calculate the letter count upto idx (not included)
    if pos < idx:
        count += k
    else:
        count -= k
    
    return count

def main():
    data=sys.argv[1]
    bw=BurrowsWheeler()
    bwt=bw.transform(data)
    inverseBwt=bw.inverse(bwt)
    print 'BWT:%s' % (bwt.replace("\0",'$'))
    print 'inverse:%s' % (inverseBwt)  

if __name__ == '__main__':
    main()
