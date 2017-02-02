# exploring the RTCGA package
# https://www.r-bloggers.com/rtcga-factory-of-r-packages-quick-guide/

# source("https://bioconductor.org/biocLite.R")
# biocLite("RTCGA")
# biocLite("RTCGA.clinical")
# biocLite("RTCGA.mRNA")
# biocLite("RTCGA.miRNASeq")
# biocLite("RTCGA.methylation")
# biocLite("RTCGA.CNV")
library("RTCGA")
library("RTCGA.clinical")
library("RTCGA.mRNA")
library("RTCGA.miRNASeq")
library("RTCGA.methylation")
library("RTCGA.CNV")

diseases <- c("BRCA", "GBM", "OV")
omics <- c("mRNA", "miRNASeq", "methylation", "CNV")

AddSurvDataToOmics(diseases, omics, save = FALSE)
