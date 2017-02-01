# exploring the RTCGA package
# https://www.r-bloggers.com/rtcga-factory-of-r-packages-quick-guide/

# source("https://bioconductor.org/biocLite.R")
# biocLite("RTCGA")
# biocLite("RTCGA.clinical")
# biocLite("RTCGA.mRNA")
# biocLite("RTCGA.miRNASeq")
# biocLite("RTCGA.methylation")
# library("RTCGA")
# library("RTCGA.clinical")
# library("RTCGA.mRNA")
# library("RTCGA.miRNASeq")
# library("RTCGA.methylation")

AddSurvivalData <- function(omics, clin){
  patients_omics <- substr(rownames(omics), 1, 12)
  omics$patient_id <- patients_omics
  merge_patients <- patients_clin %in% patients_omics
  
  patients <- patients_clin[merge_patients]
  
  days_to_death <- clin$patient.days_to_death[merge_patients]
  names(days_to_death) <- patients
  #pat <- patients[1]
  #print(days_to_death[pat])
  days_to_last_followup <- clin$patient.days_to_last_followup[merge_patients]
  names(days_to_last_followup) <- patients
  #print(days_to_last_followup)
  
  omics$days_to_death <- rep(-1, nrow(omics))
  omics$days_to_last_followup <- rep(-1, nrow(omics))

  for (patient in patients){
    omics$days_to_death[omics$patient_id == patient] <- days_to_death[patient]
    omics$days_to_last_followup[omics$patient_id == patient] <- days_to_last_followup[patient]
  }
  return(omics)
}

writeFile <- function(df, fn){
  write.table(df, sprintf("%s_with_survival_data.txt", fn))
}

# preparing clinical data
clin <- BRCA.clinical
rownames(clin) <- toupper(clin$patient.bcr_patient_barcode)
patients_clin <- rownames(clin)

# # small test
# lil <- mRNA[1:10, 1:4]
# lil2 <- AddSurvivalData(lil, clin)

# mRNA (micro array)
mRNA <- BRCA.mRNA
rownames(mRNA) <- mRNA$bcr_patient_barcode
mRNA2 <- AddSurvivalData(mRNA, clin)
writeFile(mRNA2, "mRNA2")

# miRNA-seq
miRNA <- BRCA.miRNASeq
miRNA2 <- AddSurvivalData(miRNA, clin)
writeFile(miRNA2, "miRNA2")

# methylation
methyl <- BRCA.methylation
rownames(methyl) <- toupper(methyl$bcr_patient_barcode)
methyl2 <- AddSurvivalData(methyl, clin)
writeFile(methyl2, "methylation")

# discretize days_to_death in early and late death
# tran3$days_to_death_discrete <- cut(tran3$days_to_death, breaks = 2, labels = c("early", "late"))