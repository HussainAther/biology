import re

# Cleans a DNA sequence

regex = re.compile(’ |\d|\n|\t’)
seq = ’’
for line in open(’pMOSBlue.txt’):
    seq += regex.sub(’’,line)
print(seq)
