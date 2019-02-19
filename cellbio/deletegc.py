import re

# Delete GC repeats

regex = re.compile("(?:GC){3,}")
seq="ATGATCGTACTGCGCGCTTCATGTGATGCGCGCGCGCAGACTATAAG"
print("Before:" + seq)
print("After:"+ regex.sub("",seq))
