source("https://bioconductor.org/biocLite.R")
biocLite("breastCancerNKI")
biocLite("Biobase")
library("breastCancerNKI")
library("Biobase")

data(nki)

trans_nki <- data.frame(t(exprs(nki)))
clin_nki <- pData(nki)

sum(rownames(trans_nki) != rownames(clin_nki))
# 0

# add the overall survival days to the gene expression data
trans_nki$os <- clin_nki$t.os
