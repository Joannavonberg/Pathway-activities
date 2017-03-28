diseases <- c("BRCA", "GBM", "OV", "PRAD")

# miRNAseq
for(d in diseases){
  obj <- read.table(
    sprintf(
      "gdac.broadinstitute.org_%s.Merge_mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.Level_3.2016012800.0.0/%s.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.data.txt",
      d, d),
    skip = 2)
  # assign(
  #   sprintf("%s_mRNA", d),
  #   read.table(
  #     sprintf(
  #       "/home/bit/berg0/bin/gdac.broadinstitute.org_%s.mRNA_Preprocess_Median.Level_3.2016012800.0.0/%s.medianexp.txt",
  #       d, d),
  #     skip = 2)
  # )
  con <- file(sprintf("gdac.broadinstitute.org_%s.Merge_mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.Level_3.2016012800.0.0/%s.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.data.txt", d, d))
  line <- readLines(con, 1)
  close(con)
  # obj <- get(sprintf("%s_mRNA", d))
  colnames(obj) <- strsplit(line, "\t")[[1]]
  rownames(obj) <- obj[,1]
  obj <- obj[,-1]
  assign(sprintf("%s_miRNAseq", d), data.matrix(obj))
}

# mRNA
for(d in diseases){
  obj <- read.table(
    sprintf(
      "/home/bit/berg0/bin/gdac.broadinstitute.org_%s.mRNA_Preprocess_Median.Level_3.2016012800.0.0/%s.medianexp.txt",
      d, d),
    skip = 2)
  # assign(
  #   sprintf("%s_mRNA", d),
  #   read.table(
  #     sprintf(
  #       "/home/bit/berg0/bin/gdac.broadinstitute.org_%s.mRNA_Preprocess_Median.Level_3.2016012800.0.0/%s.medianexp.txt",
  #       d, d),
  #     skip = 2)
  # )
  con <- file(sprintf("/home/bit/berg0/bin/gdac.broadinstitute.org_%s.mRNA_Preprocess_Median.Level_3.2016012800.0.0/%s.medianexp.txt", d, d))
  line <- readLines(con, 1)
  close(con)
  # obj <- get(sprintf("%s_mRNA", d))
  colnames(obj) <- strsplit(line, "\t")[[1]]
  rownames(obj) <- obj[,1]
  obj <- obj[,-1]
  assign(sprintf("%s_mRNA", d), data.matrix(obj))
}

# methylation
for(d in diseases){
  obj <- read.table(
    sprintf(
      "/home/bit/berg0/bin/gdac.broadinstitute.org_%s.Methylation_Preprocess.Level_3.2016012800.0.0/%s.meth.by_mean.data.txt",
      d, d), 
    header = TRUE,
    sep = "\t")
  # objname <- sprintf("%s_methylation", d)
  # assign(
  #   objname,
  #   read.table(
  #     sprintf(
  #       "/home/bit/berg0/bin/gdac.broadinstitute.org_%s.Methylation_Preprocess.Level_3.2016012800.0.0/%s.meth.by_mean.data.txt",
  #       d, d), 
  #     header = TRUE,
  #     sep = "\t")
  # )
  # obj <- get(objname)
  rownames(obj) <- obj[,1]
  obj[,1] <- NULL
  obj <- obj[-1,]
  assign(sprintf("%s_methylation", d), data.matrix(obj))
}

# miRNA
for(d in diseases[2]){
  obj <- read.table(
    sprintf(
      "/home/bit/berg0/Downloads/gdac.broadinstitute.org_%s.Merge_mirna__h_mirna_8x15k__unc_edu__Level_3__unc_DWD_Batch_adjusted__data.Level_3.2016012800.0.0/%s.mirna__h_mirna_8x15kv2__unc_edu__Level_3__unc_DWD_Batch_adjusted__data.data.txt",
      d, d),
    header = TRUE,
    sep = "\t")
  # objname <- sprintf("%s_methylation", d)
  # assign(
  #   objname,
  #   read.table(
  #     sprintf(
  #       "/home/bit/berg0/bin/gdac.broadinstitute.org_%s.Methylation_Preprocess.Level_3.2016012800.0.0/%s.meth.by_mean.data.txt",
  #       d, d), 
  #     header = TRUE,
  #     sep = "\t")
  # )
  # obj <- get(objname)
  rownames(obj) <- obj[,1]
  obj[,1] <- NULL
  obj <- obj[-1,]
  assign(sprintf("%s_miRNA", d), data.matrix(obj))
}

# mRNAseq
for(d in diseases){
  obj <- read.table(
    sprintf(
      "/home/bit/berg0/Downloads/gdac.broadinstitute.org_%s.mRNAseq_Preprocess.Level_3.2016012800.0.0/%s.uncv2.mRNAseq_RSEM_Z_Score.txt",
      d, d),
    header = TRUE,
    sep = "\t")
  # objname <- sprintf("%s_methylation", d)
  # assign(
  #   objname,
  #   read.table(
  #     sprintf(
  #       "/home/bit/berg0/bin/gdac.broadinstitute.org_%s.Methylation_Preprocess.Level_3.2016012800.0.0/%s.meth.by_mean.data.txt",
  #       d, d), 
  #     header = TRUE,
  #     sep = "\t")
  # )
  # obj <- get(objname)
  rownames(obj) <- obj[,1]
  obj[,1] <- NULL
  obj <- obj[-1,]
  assign(sprintf("%s_mRNAseq", d), data.matrix(obj))
}