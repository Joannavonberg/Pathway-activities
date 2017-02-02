# exploring the RTCGA package
# https://www.r-bloggers.com/rtcga-factory-of-r-packages-quick-guide/

# source("https://bioconductor.org/biocLite.R")
# biocLite("RTCGA")
# biocLite("RTCGA.clinical")
# biocLite("RTCGA.mRNA")
# biocLite("RTCGA.miRNASeq")
# biocLite("RTCGA.methylation")
library("RTCGA")
library("RTCGA.clinical")
library("RTCGA.mRNA")
library("RTCGA.miRNASeq")
library("RTCGA.methylation")



# # mRNA (micro array)
# mRNA <- BRCA.mRNA
# rownames(mRNA) <- mRNA$bcr_patient_barcode
# mRNA2 <- AddSurvivalData(mRNA, clin)
# writeFile(mRNA2, "mRNA2")
# 
# # miRNA-seq
# miRNA <- BRCA.miRNASeq
# miRNA2 <- AddSurvivalData(miRNA, clin)
# writeFile(miRNA2, "miRNA2")
# 
# # methylation
# methyl <- BRCA.methylation
# rownames(methyl) <- toupper(methyl$bcr_patient_barcode)
# methyl2 <- AddSurvivalData(methyl, clin)
# writeFile(methyl2, "methylation")

# discretize days_to_death in early and late death
# tran3$days_to_death_discrete <- cut(tran3$days_to_death, breaks = 2, labels = c("early", "late"))