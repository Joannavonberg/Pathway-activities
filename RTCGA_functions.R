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
      if(nrow(get(sprintf("%s_patients", dis_om))) == 0){
        assign(sprintf("%s_patients", dis_om), 
               NULL,
               envir = .GlobalEnv)
      }
    }
  }
}

GetPheno <- function(dis, save){
  clin <- GetTCGAData(dis, "clinical")
  df <- data.frame(
    clin[,c("patient.days_to_death", "patient.days_to_last_followup")], 
    row.names = toupper(clin$patient.bcr_patient_barcode))
  if(save){
    fn <- sprintf("%s_clinical.txt", dis)
    print(fn)
    writeFile(df, fn)
  }
  return(df)
}

GetOmicMat <- function(dis, om, save){
  df <- GetTCGAData(dis, om)
  if(is.null(df)){return(df)}
  switch(om,
         miRNASeq = NULL,
         CNV = NULL,
         rownames(df) <- toupper(df$bcr_patient_barcode)
           )
  df$bcr_patient_barcode <- NULL
  if(save){
    fn <- sprintf("%s_%s", dis, om)
    print(fn)
    writeFile(df, fn)
  }
  if(om == "mRNA"){return(t(as.matrix(df)))}
  else{return(df)}
}

GetPatientsVec <- function(barcodes, save, fn){
  res <- data.frame(patient_id = substr(barcodes, 1, 12), batch_id = as.factor(substr(barcodes, 22, 25)), row.names = toupper(barcodes))
  # res <- matrix(c(substr(barcodes, 1, 12), substr(barcodes, 22, 25)), ncol = 2, dimnames = list(barcodes))
  if(save){
    print(fn)
    writeFile(res, fn)
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

# the function writeFile wants as input:
# - a dataframe
# - an omic type for in the filename
writeFile <- function(df, fn){
  write.table(df, sprintf("%s_with_survival_data.txt", fn))
}

FindIntersection <- function(omics, disease, save = FALSE){
  # might not be a good idea, to copy very large dataframes
  oms <- mget(paste(disease, omics, "patients", sep = "_"), envir = .GlobalEnv, ifnotfound = NA)
  print(length(oms))
  print("...")
  i <- 1
  while(is.na(oms[[i]]) | is.null(oms[[i]])){
    print(i)
    i = i + 1
  }
  print("...")
  res <- oms[[i]]$patient_id
  # print(res[1:10])
  while(i <= length(oms) && !is.null(oms[[i]]) && !is.na(oms[[i]])){
    res <- res[res %in% oms[[i]]$patient_id]
    i = i + 1
    print(i)
  }
  print("......")
  if(save){
    write(unique(res), sprintf("%s_%s_patients.txt", disease, paste0(omics, collapse = "-")))
  }
  return(unique(res))
}

CorrectBatch <- function(mat, patients, pheno, save = FALSE){
  batch <- patients$batch_id
  input <- data.frame(cbind(patients, pheno[as.character(patients$patient_id),]))
  modcombat <- model.matrix(~1, data = input)
  corrected <- ComBat(dat=mat, batch=batch, mod=modcombat)
  return(corrected)
}