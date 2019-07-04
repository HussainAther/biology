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
