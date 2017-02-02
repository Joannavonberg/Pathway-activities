# This function wants as input:
# a vector of strings that indicates the which diseases should be considered (look up specific string)
# a vector of strings that indicate which omic source should be considered (look up specific string)
# a Boolean that indicates whether .txt files should be written for each omic source 
#
# and returns nothing, but creates a data-frame with survival data for each disease and omic data source combination
AddSurvDataToOmics <- function(diseases, omics, save = FALSE){
  for(dis in diseases){
    # preparing clinical data
    if(!exists("clin")){#}, where = .GlobalEnv)){
      clin <<- get(sprintf("%s.%s", dis, "clinical"))
      patients_clin <<- toupper(clin$patient.bcr_patient_barcode)
      rownames(clin) <<- patients_clin
    }
    for(om in omics){
      string <- sprintf("%s.%s", dis, om)
      if(exists(string)){
        df <- get(string)
      }
      else{
        next
      }
      print(sprintf("dimensions of df for %s: %s", om, dim(df)))
      dis_om <- sprintf("%s_%s", dis, om)
      if(om == "CNV"){
        assign(dis_om, AddSurvivalData(df, clin, om), envir = .GlobalEnv)
      }
      else{
        if(rownames(df)[1] == 1){
          rownames(df) <- toupper(df$bcr_patient_barcode)
        }
        assign(dis_om, AddSurvivalData(df, clin, om), envir = .GlobalEnv)
      }
      if(save){
        writeFile(get(dis_om), dis_om)
      }
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
AddSurvivalData <- function(omics, clin, om){
  if(om == "CNV"){
    patients_omics <- substr(omics$Sample, 1, 12)
  }
  else{
    patients_omics <- substr(rownames(omics), 1, 12)
  }
  omics$patient_id <- patients_omics
  merge_patients <- rownames(clin) %in% patients_omics
  
  patients <- rownames(clin)[merge_patients]
  days_to_death <- clin$patient.days_to_death[merge_patients]
  names(days_to_death) <- patients
  days_to_last_followup <- clin$patient.days_to_last_followup[merge_patients]
  names(days_to_last_followup) <- patients

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