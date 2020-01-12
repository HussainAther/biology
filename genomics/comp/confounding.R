library(Biobase)
library(genefilter)
library(dagdata)
data(GSE5859)

# Show an example of confounding by comparing caucations to asians
# in college admissions

eth<-factor(pData(e)$ethnicity=="CEU")
tt<-rowttests(exprs(e),eth)
HLIM<-c(0,6500)
mypar(1,2)
hist(tt$p.value,main="",xlab="p-values",nc=20,ylim=HLIM)
plot(tt$dm,-log10(tt$p.value),xlab="Effect size",ylab="-log10 (p-value)",xlim=c(-2,2))
