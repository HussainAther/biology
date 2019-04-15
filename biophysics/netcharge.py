#!/usr/bin/env python
import cgi, cgitb

"""
Web version of function to calculate the net charge of a protein and proportion of charged amino acid.
"""

def chargeandprop(AAseq):
    protseq = AAseq.upper()
    charge = -0.002
    cp = 0
    AACharge = {"C":-.045,"D":-.999,"E":-.998,"H":.091,
            "K":1,"R":1,"Y":-.001}
    for aa in protseq:
        charge += AACharge.get(aa,0)
        if aa in AACharge:
            cp += 1
        prop = float(cp)/len(AAseq)*100
    return (charge, prop)

cgitb.enable()
print("Content-Type: text/html\n")

form = cgi.FieldStorage()
uname = form.getvalue("username","NN")
seq = form.getvalue("seq","QWERTYYTREWQRTYEYTRQWE")
prop = form.getvalue("prop","n")
jobtitle = form.getvalue("title","No title")
charge, propvalue = chargeandprop(seq)
print("<html><body>Job title:"+jobtitle+"<br/>")
print("Your sequence is:<br/>"+seq+"<br/>")
print("Net charge:",charge,"<br/>")

if prop=="y":
    print("Proportion of charged AA: %.2f <br/>" %propvalue)
print "</body></html>"
