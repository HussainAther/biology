source("http://bioconductor.org/biocLite.R")
biocLite("seqLogo")
 
library(seqLogo)

"We use a probability weight matrix (PWM also called profile matrix) to describe 
the frequency of each nucleotide at each location in a motif. For the initialiation
step, we compute the frequency of each base in each position of the suspected motif. 
In the expectation step, we generatea a vector Zij that has the probability of the
motif starting in position j in sequence i. In expectation-maximization, the Z vector
gives us a way of classifying all of the nucleotides in the sequences and tells us if 
they're part of the motif. Using Bayes' rule, Zij is

Zij = Pr(Xi|Zij)Pr(Zij=1) / summation from k to L-W+1(Pr(Xi|Zij=1)Pr(Zik=1))

in which Pr(Xi|Zij=1)=Pr(Xi|Zij=1,p) is the probability of a sequence i given that
the motif starts at position j. In the maximization step, we use the Z vector to update
the PWM and the background probability distribution using 

p_c,k = (n_c,k + d_c,k) / summation (n_b,k + d_b,k)

with n_c,k = summation of i and summation of j (Zij) for k > 0 
             n_c - summation of j for n_c, j for k = 0

We repeat the expectation and maximization steps until the probability weight matrix converges. "

"Sample function for the position weight matrix."
pwm <- function(freq, total, bg=0.25){
  #using the formulae above
  p <- (freq + (sqrt(total) * 1/4)) / (total + (4 * (sqrt(total) * 1/4)))
  log2(p/bg)
}

"Calculate out all possible PWM values."
for (i in 0:8) {
  print(pwm(i, 8))
}

"Define the frequencies of nucleotides."
a <- c(0, 4, 4, 0, 3, 7, 4, 3, 5, 4, 2, 0, 0, 4)
c <- c(3, 0, 4, 8, 0, 0, 0, 3, 0, 0, 0, 0, 2, 4)
g <- c(2, 3, 0, 0, 0, 0, 0, 0, 1, 0, 6, 8, 5, 0)
t <- c(3, 1, 0, 0, 5, 1, 4, 2, 2, 4, 0, 0, 1, 0)
m <- matrix(data=c(a,c,g,t), nrow=4, byrow=T, dimnames=list(c("a", "c", "g", "t")))

"Add the matrix to our pwm."
mm <- pwm(m, 8)

"For some sample sequence."
seq <- "ttacataagtagtc"
x <- strsplit(x=seq, split="")

"Initialize a vector."
seq_score <- vector()

"Get the corresponding values."
for (i in 1:14) {
  seq_score[i] <- mm[x[[1]][i],i]
}

"Retrieve the max score."
sum(apply(mm, 2, max))

