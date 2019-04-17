import re

# Delete GC repeats using re

regex = re.compile("(?:GC){3,}")
seq="ATGATCGTACTGCGCGCTTCATGTGATGCGCGCGCGCAGACTATAAG"
print("Before:" + seq)
print("After:"+ regex.sub("",seq))
