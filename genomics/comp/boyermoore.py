"""
The Boyer-Moore (boyer moore Boyer Moore) uses bad character heuristic and good suffix
heuristic to search for patterns within strings. Both of the above heuristics can also 
be used independently to search a pattern in a text. Let us first understand how two 
independent approaches work together in the Boyer Moore algorithm. If we take a look 
at the Naive algorithm, it slides the pattern over the text one by one. KMP algorithm 
does preprocessing over the pattern so that the pattern can be shifted by more than one. 
The Boyer Moore algorithm does preprocessing for the same reason. It preporcesses the 
pattern and creates different arrays for both heuristics. At every step, it slides the 
pattern by max of the slides suggested by the two heuristics. So it uses best of the 
two heuristics at every step.
"""

def badCharHeuristic(string, size): 
    """ 
    The preprocessing function for 
    Boyer Moore's bad character heuristic 
    """
    badChar = [-1] * 256 # Initialize all occurence as -1 
    for i in range(size): # Fill the actual value of last occurence 
        badChar[ord(string[i])] = i; 
    return badChar # return the list  
  
def search(txt, pat): 
    """ 
    A pattern searching function that uses Bad Character 
    Heuristic of Boyer Moore Algorithm 
    """
    m = len(pat) 
    n = len(txt) 
    # create the bad character list by calling  
    # the preprocessing function badCharHeuristic() 
    # for given pattern 
    badChar = badCharHeuristic(pat, m)  
    # s is shift of the pattern with respect to text 
    s = 0
    while(s <= n-m): 
        j = m-1
        # Keep reducing index j of pattern while  
        # characters of pattern and text are matching 
        # at this shift s 
        while j>=0 and pat[j] == txt[s+j]: 
            j -= 1
        # If the pattern is present at current shift,  
        # then index j will become -1 after the above loop 
        if j<0: 
            print("Pattern occur at shift = {}".format(s)) 
            """ 
            Shift the pattern so that the next character in text 
            aligns with the last occurrence of it in pattern. 
            The condition s+m < n is necessary for the case when 
            pattern occurs at the end of text 
            """
            s += (m-badChar[ord(txt[s+m])] if s+m<n else 1) 
        else: 
            """ 
            Shift the pattern so that the bad character in text 
            aligns with the last occurrence of it in pattern. The 
            max function is used to make sure that we get a positive 
            shift. We may get a negative shift if the last occurrence 
            of bad character in pattern is on the right side of the 
            current character. 
            """
            s += max(1, j-badChar[ord(txt[s+j])]) 
