# exploring the RTCGA package
# https://www.r-bloggers.com/rtcga-factory-of-r-packages-quick-guide/

source("https://bioconductor.org/biocLite.R")
#biocLite("RTCGA")
biocLite("RTCGA.clinical")
biocLite("RTCGA.mRNA")
library("RTCGA")
library("RTCGA.clinical")
library("RTCGA.mRNA")
library("ggplot2")
library("survminer")

options(stringsAsFactors = FALSE)

names <- checkTCGA(what = "DataSets", cancerType = "BRCA")

# downloading clinical data and making a Kaplan-Meier plot
#dataset <- downloadTCGA("BRCA", dataSet = "Merge_Clinical.Level_1", destDir = ".")
clinical_data <- readTCGA(path = file.path("./gdac.broadinstitute.org_BRCA.Merge_Clinical.Level_1.2016012800.0.0/BRCA.clin.merged.txt"), dataType = "clinical")

survival <- survivalTCGA(data)
kmTCGA(survival)

# downloading transcriptomic data
#dataset <- downloadTCGA("BRCA", dataSet = names[34,"Name"], destDir = ".")
#path <- file.path("gdac.broadinstitute.org_BRCA.Merge_transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.Level_3.2016012800.0.0/BRCA.transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.data")
path <- file.path("./gdac.broadinstitute.org_BRCA.Merge_transcriptome__agilent.2016012800.0.0/BRCA.transcriptome__agilent.data.txt")
transcriptomic_data <- readTCGA(path = path, dataType = "mRNA")

# connecting the two databases, using barcode as identifier (since the UUID is only available for the clinical data)
# https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode?preview=/39294833/39421313/creating_barcodes.png
clin2 <- clinical_data
rownames(clin2) <- toupper(clin2$patient.bcr_patient_barcode)
tran2 <- transcriptomic_data
rownames(tran2) <- toupper(tran2$bcr_patient_barcode)

# patients_clin <- rownames(clin2)
# patients_tran <- unique(substr(rownames(tran2), 1, 12))
# merge_patients <- patients_clin %in% patients_tran
# # merge_patients <- patients_tran %in% patients_clin

patients_clin <- rownames(clin2)
patients_tran <- substr(rownames(tran2), 1, 12)
tran2$patient_id <- patients_tran
merge_patients <- patients_clin %in% patients_tran
#merge_patients <- patients_tran %in% patients_clin

# keep only the patients with clinical data
nrow(tran2)
# 590
sum(merge_patients)
# 527
# there's a difference, because there are more samples than patients (multiple samples are possible per patient)
days_to_death <- clin2$patient.days_to_death[merge_patients]
names(days_to_death) <- rownames(clin2)[merge_patients]

# add the days to death for each transcriptomic sample
tran2$days_to_death <- rep(-1, nrow(tran2))
for (patient in patients_clin){
  tran2$days_to_death[tran2$patient_id == patient] <- days_to_death[patient]
}

# remove all samples for which days_to_death == NA (better to do this before, if it is at all sensible to do)
tran3 <- tran2[!is.na(tran2$days_to_death),]

# discretize days_to_death in early and late death
tran3$days_to_death_discrete <- cut(tran3$days_to_death, breaks = 2, labels = c("early", "late"))

# integrate gene lists based on pathways to perform feature selection

# download pathway gene lists