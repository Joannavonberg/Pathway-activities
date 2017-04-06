diseases <- c("BRCA", "GBM", "OV", "PRAD")
omics <- c("mRNA", "mRNASeq", "CNV", "miRNA", "miRNASeq", "methylation")

for(d in diseases){
  load(sprintf("/home/abidata/Joanna/%s/%s_pheno.rda", d, d))
  for(o in omics){
    try(
      load(sprintf("/home/abidata/Joanna/%s/%s/%s_%s.rda", d, o, d, o))
    )
    try(
      load(sprintf("/home/abidata/Joanna/%s/%s/%s_%s_patients.rda", d, o, d, o))
    )
  }
}