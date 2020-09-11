#' perform GWAS
#' @param phenos a vector
#'@export

perform_gwas <- function(phenos, sim_label, geno_path){
  # create phenetype file that has a format suitable for plink
  prefix <- paste0("output/sim", sim_label, "/")
  val_pheno_file <- paste0(prefix, "val.pheno")
  fam_pheno_file <- paste0(prefix, "fam.pheno")
  pheno_file <- paste0(prefix, "pheno.pheno")
  gwas_res_file <- paste0(prefix, "gwas_res")
  write.table(phenos, file=val_pheno_file ,row.names=FALSE, col.names = FALSE)
  fam_file <- paste0(geno_path, ".fam")
  system(paste("awk '{print $1,$2}'", fam_file, ">", fam_pheno_file))
  system(paste("paste", fam_pheno_file, val_pheno_file, ">", pheno_file))
  # perform plink association tests to get p values
  system(paste("plink --bfile",  geno_path, "--assoc --pheno", pheno_file, "--allow-no-sex --out", gwas_res_file))
  # gwas_res is a df with SNP names and p values
  gwas_df <- read.table("output/sim1/gwas_res.qassoc", head=T)
  return(gwas_df)
}
