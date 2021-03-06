# This function wants as input:
# a vector of strings that indicates the which diseases should be considered (look up specific string)
# a vector of strings that indicate which omic source should be considered (look up specific string)
# a Boolean that indicates whether .txt files should be written for each omic source 
#
# and returns nothing, but creates a data-frame with survival data for each disease and omic data source combination
AddSurvDataToOmics <- function(diseases, omics, save = FALSE){
  for(dis in names(diseases)){
    res <- GetPheno(dis, save, diseases[dis])
    assign(sprintf("%s_pheno", dis), 
           res[[1]],
           envir = .GlobalEnv)
    for(om in omics){
      dis_om <- sprintf("%s_%s", dis, om)
      assign(dis_om, 
             GetOmicMat(dis, om, save, diseases[dis]), 
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

GetPheno <- function(dis, save, newdata){
  res <- GetTCGAData(dis, "clinical", newdata)
  clin <- res[[1]]
  # df <- data.frame(
  #   clin[,c("patient.days_to_death", "patient.days_to_last_followup")], 
  #   row.names = toupper(clin$patient.bcr_patient_barcode))
  days <- ifelse(is.na(clin$patient.days_to_death), 
                 clin$patient.days_to_last_followup, 
                 clin$patient.days_to_death)
  censored <- ifelse(is.na(clin$patient.days_to_death),
                     TRUE,
                     FALSE)
  df <- data.frame(cbind(days, censored), row.names = clin$patient.bcr_patient_barcode)
  if(save){
    fn <- sprintf("%s_clinical.txt", dis)
    writeFile(df, fn)
  }
  return(list(df, res[[2]]))
}

NewPheno <- function(old_pheno){
  days <- ifelse(is.na(old_pheno$patient.days_to_death), 
                 old_pheno$patient.days_to_last_followup, 
                 old_pheno$patient.days_to_death)
  censored <- ifelse(is.na(old_pheno$patient.days_to_death),
                     TRUE,
                     FALSE)

  return()
}

GetOmicMat <- function(dis, om, save, newdata){
  df <- GetTCGAData(dis, om, newdata)[[1]]
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
  res <- data.frame(patient_id = substr(barcodes, 1, 12), batch_id = as.factor(substr(barcodes, 6, 7)), row.names = toupper(barcodes))
  # res <- matrix(c(substr(barcodes, 1, 12), substr(barcodes, 22, 25)), ncol = 2, dimnames = list(barcodes))
  if(save){
    print(fn)
    writeFile(res, fn)
  }
  return(res)
}

GetTCGAData <- function(dis, type, newdata = TRUE){
  if(newdata){
    res <- sprintf("%s.%s.20160128", dis, type)
    if(exists(res)){
      return(list(get(res), TRUE))
    }
  }
  else{
    res <- sprintf("%s.%s", dis, type)
    if(exists(res)){
      return(list(get(res), FALSE))
    }
    return(NULL)
  }
  return(NULL)
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
  # print(length(oms))
  # print("...")
  i <- 1
  while(is.null(oms[[i]]) || is.na(oms[[i]])){
    # print(i)
    i = i + 1
  }
  # print("...")
  res <- as.character(oms[[i]]$patient_id)
  # print(res[1:4])
  # print(res[1:10])
  while(i <= length(oms) && !is.null(oms[[i]]) && !is.na(oms[[i]])){
    res <- res[res %in% as.character(oms[[i]]$patient_id)]
    # print(res[1:2])
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
  ind <- rownames(pheno) %in% as.character(patients$patient_id)
  # cbind(patients, pheno[ind, "days"])
  res <- patients
  res$pheno <- rep(-1, nrow(res))
  for(i in rownames(pheno)[ind]){
    res[res$patient_id == i, "pheno"] <- pheno[i, "days"]
  }
  input <- data.frame(res)
  modcombat <- model.matrix(~1, data = input)
  corrected <- ComBat(dat=mat, batch=batch, mod=modcombat)
  return(corrected)
}

MakeBatchPlot <- function(){
  # pc <- prcomp(t(BRCA_mRNA_imputed$data))
  # pc.pred <- predict(pc,newdata = t(BRCA_mRNA_imputed$data))
  p <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour = as.numeric(new_pheno[BRCA_mRNA_patients$patient_id,"days"]))) + ggtitle("Batch Effects before Transformation") + labs(y = "PC1", x = "PC2", colour = "Classes") + geom_point(shape=19)
  p
  return(p)
}

PathwaysMrna <- Vectorize(function(i, dis, fun = mean()){
  obj <- get(sprintf("%s_mRNA", dis))
  ind <- rownames(obj) %in% names(entrez_mRNA[as.numeric(kegg2[[i]])])
  if(sum(ind) == 0){return(NULL)}
  if(sum(ind) == 1){return(obj[ind,])}
  return(apply(obj[ind,], 2, fun))
})

FirstPCPerPathway <- function(om, data, patients, lookup = entrez_hgnc, pwlist = kegg2, mseq = FALSE){
  if(mseq){
    genes <- rownames(data)
  }
  if(om == "miRNA" | om == "miRNASeq"){
    genes <- rownames(data)
  }
  else{
    genes <- lookup[rownames(data), "entrezgene"]
    pwlist <- pwlist[[1]]
  }
  res <- do.call(rbind, 
                 lapply(
                   pwlist,
                   pc1_per_pathway,
                   genes,
                   om,
                   data
                 )
  )
  colnames(res) <- colnames(data)
  return(res)
}

pc1_per_pathway <- function(pw_genes, genes, om, data){
  ind <- genes %in% pw_genes
  if(om == "CNV"){fun <- abs}else{fun <- identity}
  if(sum(ind) == 0){return(NULL)}
  if(sum(ind) == 1){return(fun(data[ind,]))}
  ind2 <- fun(data[ind,])
  pw2 <- predict(prcomp(t(ind2), scale. = TRUE), newdata = t(ind2))[,1]
  names(pw2) <- NULL
  return(pw2)
}

HgncToEntrez <- function(hgnc){
  return(entrez_hgnc[hgnc,'entrezgene'])
}

EntrezToHgnc <- function(entrez){
  return(entrez_hgnc[entrez_hgnc$entrezgene %in% entrez, 'hgnc_symbol'])
}

RUCN <- function(data, patients){
  genes <- HgncToEntrez(rownames(data))
  res <- do.call(rbind, 
                 lapply(
                   kegg2,
                   function(pw_genes){
                     overlap_entrez <- genes[genes %in% pw_genes]
                     if(length(overlap_entrez) == 0){return(NULL)}
                     overlap_hgnc <- EntrezToHgnc(overlap_entrez)
                     pw_data <- data[overlap_hgnc,]
                     if(length(overlap_entrez) == 1){return(abs(pw_data)/length(pw_genes))}
                     nr <- nrow(pw_data)
                     pw_ratio <- apply(
                       pw_data, 2, function(CNV_pw){
                         return(mean(abs(CNV_pw)))
                       }
                     )
                     names(pw_ratio) <- NULL
                     return(pw_ratio)
                   }
                 )
  )
  colnames(res) <- rownames(patients)
  return(res)
}

avgperpw <- function(i, dis, om, fun = mean, top50 = FALSE, pval = NULL){
  # print("i is:")
  # print(i)
  # print("\n")
  obj <- get(sprintf("%s_%s", dis, om))
  fun2 <- identity
  fun3 <- fun2
  if(om == "CNV"){
    fun2 <- abs
  }
  if(om == "miRNA" | om == "miRNAseq"){
    ind <- rownames(obj) %in% pw_miRNA[[i]]
    fun3 <- HgncToEntrez
  }
  else{
    ind <- rownames(obj) %in% EntrezToHgnc(as.numeric(kegg2[[i]]))
    # print(sum(ind))
  }
  if(sum(ind) == 0){return(NULL)}
  if(sum(ind) == 1){return(fun2(obj[ind,]))}
  if(top50){
    # print(pval[1:4])
    pval_pw <- pval[rownames(obj)[ind]]
    # print(rownames(obj)[ind])
    # print(pval_pw)
    order1 <- order(pval_pw, decreasing = FALSE)
    # print(order1)
    top <- names(pval_pw[order1[1:ceiling(length(order1)/2)]])
    # print(top)
    if(length(top) == 1){return(fun2(obj[fun3(top),]))}
    return(apply(fun2(obj[fun3(top),]), 2, fun))
  }
  else{
    return(apply(fun2(obj[ind,]), 2, fun))
  }
}

AveragePerPathway <- Vectorize(avgperpw, vectorize.args = "i")
