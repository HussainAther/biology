"""
Read enzyme data.
"""

def read_enzymes_from_file(filename):
    """
    Read enzymes.
    """
    with open(filename) as file:
        skip_intro(file)
        return (get_enzymes(file), get_references(file))

def skip_intro(file):
    """
    Skip through the documentation that appears at the beginning
    of file, leaving it positioned at the first line of the first enzyme.
    """
    line = ""
    while not line.startswith('<REFERENCES>'):
        line = file.readline()
    while len(line) > 1: # always 1 for '\n'
        line = file.readline()
    return line
