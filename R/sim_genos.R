#' Generate genotypes for three populations YRI, CEU, and CHB.
#' Note: to run this function, make sure system or terminal or cluster has commands stdpopsim, tskit, and plink2
#' usually on the PC, make sure these commands have to be under /usr/local/bin
#' 
#' @param N_YRI the number of YRI individuals
#' @param N_CEU the number of CEU individuals
#' @param N_CHB the number of CHB individuals
#' @param sim_label denotes the x-th simulation
#' @return output: path for genotypes, a string
#' @return freq_stat: path for MAF frequency statistic summary
#' @export

sim_genos <- function(N_YRI, N_CEU, N_CHB, sim_label=1){
  # create a directory for simulation output
  prefix <- paste0("output/sim", sim_label, "/")
  system(paste("mkdir -p", prefix))
  output <- paste0(prefix, "chr20Y", N_YRI, "E", N_CEU, "A", N_CHB)
  output_ts <- paste0(output, ".ts")
  output_vcf <- paste0(output, ".vcf")
  # simulate using out-of-Africa model wrapped by stdpopsim
  system(paste("stdpopsim HomSap -s 1046 -g HapMapII_GRCh37 -c chr20 -o", output_ts, "-d OutOfAfrica_3G09", N_YRI, N_CEU, N_CHB))
  # convert output to a vcf file
  system(paste("tskit vcf", output_ts, ">", output_vcf))
  # convert to a bed file
  system(paste("plink2 -vcf", output_vcf, "--make-bed  --out", output))
  # use plink2 for QC
  # remove snps with maf < 0.01
  system(paste("plink2 --bfile", output, "--maf 0.01 --make-bed --out", output))
  # remove genotyping error
  system(paste("plink2 --bfile", output, "--hwe 1e-50 keep-fewhet --make-bed --out", output))
  # remove variants that fail HWE but still keep the population stratification
  system(paste("plink2 --bfile", output, "--hwe 1e-5 keep-fewhet --make-bed --out", output))
  # also generate statistic of MAF file
  freq_stat <- paste0(prefix,"freq_stat")
  system(paste("plink2 --bfile", output, "--freq --out", freq_stat))
  
  return(list(output, freq_stat)) # note output and freq_stat are strings
}

#res<- sim_genos(10, 10, 10)