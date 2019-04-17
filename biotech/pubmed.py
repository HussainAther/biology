from Bio import Entrez

"""
Retrieve and display data from PubMed.
"""

my_em = "user@example.com"
db = "gene"
term = "cobalamin synthase homo sapiens"
h_search = Entrez.esearch(db=db, email=my_em, term=term)
record = Entrez.read(h_search)
res_ids = record["IdList"]
for r_id in res_ids:
    h_summ = Entrez.esummary(db=db, id=r_id, email=my_em)
    summ = Entrez.read(h_summ)
    print(r_id)
    print(summ[0]["Description"])
    print(summ[0]["Summary"])
    print("==============================================")
