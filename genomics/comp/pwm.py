"""
We use a probability weight matrix (PWM also called profile matrix) to describe 
the frequency of each nucleotide at each location in a motif. For the initialiation
step, we compute the frequency of each base in each position of the suspected motif. 
In the expectation step, we generatea a vector Zij that has the probability of the
motif starting in position j in sequence i. In expectation-maximization, the Z vector
gives us a way of classifying all of the nucleotides in the sequences and tells us if 
they're part of the motif. Using Bayes' rule, Zij is

Zij = Pr(Xi|Zij)Pr(Zij=1) / summation from k to L-W+1(Pr(Xi|Zij=1)Pr(Zik=1))

in which Pr(Xi|Zij=1)=Pr(Xi|Zij=1,p) is the probability of a sequence i given that
the motif starts at position j. 
"""
