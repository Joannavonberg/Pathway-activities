# This function wants as input:
# a string that indicates the disease (look up specific string)
# a vector of strings that indicate which omic source should be considered (look up specific string)
# a Boolean that indicates whether .txt files should be written for each omic source 
#
# and returns nothing, but it creates a data-frame with survival data for each omic data source
AddSurvDataToOmics <- function(disease, omics, save = FALSE){
  # preparing clinical data
  if(!exists("clin")){
    clin <- get(sprintf("%s.%s", disease, "clinical"))
    rownames(clin) <- toupper(clin$patient.bcr_patient_barcode)
    patients_clin <- rownames(clin)
  }
  for(om in omics){
    assign(om, get(sprintf("%s.%s", disease, om)))
    df <- get(om)
    rownames(df) <- toupper(df$bcr_patient_barcode)
    df <- AddSurvivalData(df, clin)
    if(save){
      writeFile(df, om)
    }
  }
}

# the function AddSurvivalData wants as input: 
# - a dataframe with samples as rows and omic features of one omic source as columns (omics)
# - a dataframe with patients as rows and clinical features as columns
#
# and returns:
# - a dataframe like the original omic input, with extra columns patient_id, days_to_death and days_to_last_followup
# NB. It is not specific to any disease
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

# the function writeFile wants as input:
# - a dataframe
# - an omic type for in the filename
writeFile <- function(df, fn){
  write.table(df, sprintf("%s_with_survival_data.txt", fn))
}