library(knitr)
library(GenomicFeatures)
library(GenomicRanges)
library(IRanges)
library(devtools)
library(ERBS)
library(Homo.sapiens)
opts_chunk$set(fig.path=paste0("figure/", sub("(.*).Rmd","\\1",basename(knitr:::knit_concord$get('infile'))), "-"))

# Genomic sequence utility for motif tracking
