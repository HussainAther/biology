from Bio import Entrez
from Bio import Medline
from tqdm import tqdm

"""
Biopython work for fetching abstracts.
"""

# Change this email to your email address
Entrez.email = "shussainather@gmail.com"

keyword = "optical trap"

result = Entrez.read(Entrez.esearch(db="pubmed", retmax=10, term=keyword))
print(
    "Total number of publications that contain the term {}: {}".format(
        keyword, result["Count"]
    )
)

