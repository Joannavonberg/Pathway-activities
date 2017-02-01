# using pooled data from 8 data sets taken together
# see van Vliet et al, 2008

# install.packages("R.matlab")
library("R.matlab")

data <- readMat("van_Vliet_2008_data/DMMPLCN_080423.mat")

SampleLabels <- unlist(data[[1]][[4]])
ReporterIDs <- unlist(data[[1]][[5]])

test <- unlist(data[[1]][[3]])
test2 <- as.data.frame(t(test))
rownames(test2) <- SampleLabels
colnames(test2) <- ReporterIDs

genex <- test2


clin <- unlist(data[[1]][[7]])
labs <- unlist(clin[1])
