library(knitr)
library(GenomicFeatures)
library(GenomicRanges)
library(IRanges)
library(devtools)
library(ERBS)
library(Homo.sapiens)
opts_chunk$set(fig.path=paste0("figure/", sub("(.*).Rmd","\\1",basename(knitr:::knit_concord$get('infile'))), "-"))

# Genomic sequence utility for motif tracking

library(ERBS)
data(HepG2)
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens

# Reference sequence
Hsapiens$chr17

# Targeted retrieval of reference sequence
hepseq = getSeq(Hsapiens, HepG2)
rhepseq = getSeq(Hsapiens, shift(HepG2,2500))
hepseq

# Count motif occurrences 
sum(vcountPattern("TCAAGGTCA", hepseq))+sum(vcountPattern("TCAAGGTCA", 
   reverseComplement(hepseq)))
sum(vcountPattern("TCAAGGTCA", rhepseq))+sum(vcountPattern("TCAAGGTCA", 
   reverseComplement(rhepseq)))
