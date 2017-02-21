# This function wants as input:
# a vector of strings that indicates the which diseases should be considered (look up specific string)
# a vector of strings that indicate which omic source should be considered (look up specific string)
# a Boolean that indicates whether .txt files should be written for each omic source 
#
# and returns nothing, but creates a data-frame with survival data for each disease and omic data source combination
AddSurvDataToOmics <- function(diseases, omics, save = FALSE){
  for(dis in diseases){
    assign(sprintf("%s_pheno", dis), 
           GetPheno(dis, save),
           envir = .GlobalEnv)
    for(om in omics){
      dis_om <- sprintf("%s_%s", dis, om)
      assign(dis_om, 
             GetOmicMat(dis, om, save), 
             envir = .GlobalEnv)
      assign(sprintf("%s_patients", dis_om), 
             GetPatientsVec(colnames(get(dis_om)), save, sprintf("%s_patients.txt", dis_om)), 
             envir = .GlobalEnv)
    }
  }
}

GetPheno <- function(dis, save){
  clin <- GetTCGAData(dis, "clinical")
  df <- data.frame(
    clin[,c("patient.days_to_death", "patient.days_to_last_followup")], 
    row.names = toupper(clin$patient.bcr_patient_barcode))
  if(save){
    writeFile(df, sprintf("%s_clinical.txt", dis))
  }
  return(df)
}

GetOmicMat <- function(dis, om, save){
  df <- GetTCGAData(dis, om)
  rownames(df) <- toupper(df$bcr_patient_barcode)
  df$bcr_patient_barcode <- NULL
  if(save){
    writeFile(df, dis_om)
  }
  return(t(as.matrix(df)))
}

GetPatientsVec <- function(barcodes, save, fn){
  res <- data.frame(patient_id = substr(barcodes, 1, 12), batch_id = as.factor(substr(barcodes, 22, 25)), row.names = toupper(barcodes))
  # res <- matrix(c(substr(barcodes, 1, 12), substr(barcodes, 22, 25)), ncol = 2, dimnames = list(barcodes))
  if(save){
    write(res, fn, ncolumns = 1)
  }
  return(res)
}

GetTCGAData <- function(dis, type){
  string <- sprintf("%s.%s.20160128", dis, type)
  if(exists(string)){
    return(get(string))
  }
  else{
    string <- sprintf("%s.%s", dis, type)
    if(exists(string)){
      return(get(string))
    }
    else{
      return(NULL)
    }
  }
}

# PhenoData <- function(pheno, dis){
#     # preparing clinical data
#     new_data <- exists(sprintf("%s.%s.20160128", dis, "clinical"))
#     if(new_data){
#       clin <<- get(sprintf("%s.%s.20160128", dis, "clinical"), envir = .GlobalEnv)
#       print(sprintf("For %s, clinical data from 2016-01-28 will be used", dis))
#     }
#     else{
#       clin <<- get(sprintf("%s.%s", dis, "clinical"), envir = .GlobalEnv)
#       print(sprintf("Clinical data from 2016-01-28 is not available, will use old data for %s.", dis))
#     }
#     barcodes_clin <<- toupper(clin$patient.bcr_patient_barcode)
#     print(barcodes_clin[1:4])
#     rownames(clin) <<- barcodes_clin
#     # patients <- substr(barcodes_clin, 1, 12)
#     batch_id <- substr(barcodes_clin, 22, 25)
# 
#     days_to_death <- clin$patient.days_to_death
#     # names(days_to_death) <- 
#     days_to_last_followup <- clin$patient.days_to_last_followup
#     # names(days_to_last_followup) <- patients
#     return(cbind(barcodes_clin, batch_id, days_to_death, days_to_last_followup))
# looks like:
# barcode - batch_id - days_to_death - days_to_last_dinges
# }

# the function AddSurvivalData wants as input: 
# - a dataframe with samples as rows and omic features of one omic source as columns (omics)
# - a dataframe with patients as rows and clinical features as columns
#
# and returns:
# - a dataframe like the original omic input, with extra columns patient_id, days_to_death and days_to_last_followup
# NB. It is not specific to any disease
# AddSurvivalData <- function(omics, clin, om, prad = FALSE){
#   if(om == "CNV"){
#     patients_omics <- substr(omics$Sample, 1, 12)
#     batch_id <- substr(omics$Sample, 22, 25)
#   }
#   else{
#     patients_omics <- substr(rownames(omics), 1, 12)
#     batch_id <- substr(rownames(omics), 22, 25)
#   }
#   omics$batch_id <- batch_id
#   omics$patient_id <- patients_omics
#   merge_patients <- rownames(clin) %in% patients_omics
#   
#   patients <- rownames(clin)[merge_patients]
#   days_to_death <- clin$patient.days_to_death[merge_patients]
#   names(days_to_death) <- patients
#   days_to_last_followup <- clin$patient.days_to_last_followup[merge_patients]
#   names(days_to_last_followup) <- patients
#   
#   omics$days_to_death <- rep(-1, nrow(omics))
#   omics$days_to_last_followup <- rep(-1, nrow(omics))
#   for (patient in patients){
#     omics$days_to_death[omics$patient_id == patient] <- days_to_death[patient]
#     omics$days_to_last_followup[omics$patient_id == patient] <- days_to_last_followup[patient]
#   }
#   return(omics)
# }

# the function writeFile wants as input:
# - a dataframe
# - an omic type for in the filename
writeFile <- function(df, fn){
  write.table(df, sprintf("%s_with_survival_data.txt", fn))
}

FindIntersection <- function(omics, disease, save = FALSE){
  # might not be a good idea, to copy very large dataframes
  oms <- mget(paste(disease, omics, sep = "_"), envir = .GlobalEnv, ifnotfound = NA)
  # print(length(oms))
  # print("...")
  i <- 1
  while(is.na(oms[[i]])){
    # print(i)
    i = i + 1
  }
  # print("...")
  res <- oms[[i]]$patient_id
  # print(res[1:10])
  while(i <= length(oms) && !is.na(oms[[i]])){
    res <- res[res %in% oms[[i]]$patient_id]
    i = i + 1
    # print(i)
  }
  # print("......")
  if(save){
    write(unique(res), sprintf("%s_%s_patients.txt", disease, paste0(omics, collapse = "-")))
  }
  return(unique(res))
}

CorrectBatch <- function(mat, patients, pheno, save = FALSE){
  batch <- patients$batch_id
  modcombat <- model.matrix(~1, data = pheno)
  corrected <- ComBat(dat=mat, batch=batch, mod=modcombat)
  return(corrected)
}