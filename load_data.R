diseases <- c("BRCA", "GBM", "OV")
omics <- c("mRNA", "CNV", "miRNASeq", "methylation")

for(d in diseases){
  load(sprintf("/home/abidata/Joanna/%s/%s_pheno.rda", d, d))
  for(o in omics){
    load(sprintf("/home/abidata/Joanna/%s/%s/%s_%s.rda", d, o, d, o))
    load(sprintf("/home/abidata/Joanna/%s/%s/%s_%s_patients.rda", d, o, d, o))
  }
}