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

diseases <- c("BRCA")#, "GBM", "OV")
omics <- c("mRNA", "miRNASeq", "methylation", "CNV")

AddSurvDataToOmics(diseases, omics, save = FALSE)

for(dis in diseases){
  assign(sprintf("%s_inter", dis), FindIntersection(omics, dis, save = TRUE))
}

# problem: for GBM, there is only miRNA-seq data for 5 patients and these do not show up in both the CNV and the methylation data
# therefore there are no patients for which all data is there

GBM_miRNASeq_donotuse <- GBM_miRNASeq
GBM_miRNASeq <- NA

GBM_inter <- FindIntersection(omics, "GBM")
# there are 283 GBM patients for which there is methylation and CNV data
# I should probably look for mRNA and more miRNA data