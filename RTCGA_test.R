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

# next step: see for which patients all data types are available, and save just the patient barcode's (no need to save all the data seperately)

t1 <- 1:10
t2 <- 4:20

# logical vector:
ind1 <- t1 %in% t2
# values that are TRUE:
which(t1 %in% t2)

t3 <- c(2, 6, 18, 22)

# this doesn't work:
ind2 <- t1 %in% t2 %in% t3

# test:

# which patients are in t1 and t2?
# NB: which returns INDICES!
t1_and_t2 <- which(t1 %in% t2)
t1_and_t2_and_t3 <- which(t3 %in% t1[t1_and_t2])

# testing with real data

m_test <- BRCA_mRNA[1:20, 1:4]
m_test$patient_id <- BRCA_mRNA$patient_id[1:20]
mi_test <- BRCA_miRNASeq[1:20, 1:4]
mi_test$patient_id <- BRCA_miRNASeq$patient_id[1:20]

m_mi <- which(m_test$patient_id %in% mi_test$patient_id)

me_test <- BRCA_methylation[1:20, 1:4]
me_test$patient_id <- BRCA_methylation$patient_id[1:20]

m_mi_me <- which(me_test$patient_id %in% m_test$patient_id[m_mi])
# returns an empty vector, but I think that is correct

# testing with the whole dataframe

mrna_mirna <- which(BRCA_mRNA$patient_id %in% BRCA_miRNASeq$patient_id)
length(unique(BRCA_mRNA[mrna_mirna, "patient_id"]))
# 512
mrna_mirna_methyl <- which(BRCA_methylation$patient_id %in% BRCA_mRNA$patient_id[mrna_mirna])
length(unique(BRCA_methylation$patient_id[mrna_mirna_methyl]))
# 302
mrna_mirna_methyl_cnv <- which(BRCA_CNV$patient_id %in% BRCA_methylation$patient_id[mrna_mirna_methyl]) 
length(unique(BRCA_CNV[mrna_mirna_methyl_cnv, "patient_id"]))
# 302

length(mrna_mirna_methyl_cnv)
# 91169??

length(unique(mrna_mirna_methyl_cnv))
# still 91169??