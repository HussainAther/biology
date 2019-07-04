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
    while len(line) > 1: # always 1 for "\n"
        line = file.readline()
    return line

def get_enzymes(src):
    """
    Get the enzymes.
    """
    enzymes = {}
    enzyme = next_enzyme(src)
    while enzyme:
        enzymes[enzyme[0]] = enzyme # dict key is enzyme's name
        enzyme = next_enzyme(src)
    return enzymes

def read_field(file):
    return file.readline()[3:-1]

def next_enzyme(file):
    """
    Read the data for the next enzyme, returning a list of the
    form: [enzyme_name, prototype, source, recognition_tuple,
    (genus, species, subspecies), references_tuple].
    """
    name = read_field(file)
    if name: # otherwise last enzyme read
        fields = [name] + read_other_fields(file)
        fields[2] = parse_organism(fields[2])
        fields[7] = [int(num) for num in fields[7].split(',')]
        file.readline() # skip blank line
        return fields
