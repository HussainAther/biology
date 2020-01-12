library(knitr)
opts_chunk$set(fig.path=paste0("figure/", sub("(.*).Rmd","\\1",basename(knitr:::knit_concord$get('infile'))), "-"))

# Model batch effects with factor analysis

library(GSE5859Subset)
data(GSE5859Subset)
