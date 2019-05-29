"""
Sequence by hybridization (SBH) sequencing approach that uses a miniature DNA array
of a DNA chip that has thousands of short DNA fragments called probes. Given a short 
probe (an 8- to 30-nucleotide single-stranded synthetic DNA fragment) and a single-stranded 
target DNA fragment, the target will hybridize with the probe if the probe is a substring 
of the target’s Watson-Crick complement. When the probe and the target are mixed together, 
they form a weak chemical bond and stick together. For example, a probe ACCGTGGA
will hybridize to CCCTGGCACCTA since it's complementary to the substring TGGCACCT of the target.
"""

def sbh(x, s):
    """
    For a list x of all l-mers of a string s, return string s such that Spectrum(s, l) = x. 
    """
