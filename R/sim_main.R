#' main function for simulation
#'
#' @param sim_label denotes the x-th simulation
#' @param N_YRI the number of YRI individuals
#' @param N_CEU the number of CEU individuals
#' @param N_CHB the number of CHB individuals
#' @param p_YRI the proportion of causal SNPs in YRI population
#' @param p_CEU the proportion of causal SNPs in CEU population
#' @param p_CHB the proportion of causal SNPs in CHB population
#' @param h2_YRI heritability of this phenotype in YRI population
#' @param h2_CEU heritability of this phenotype in CEU population
#' @param h2_CHB heritability of this phenotype in CHB population
#' @param prare_YRI proportion of causal rare SNPs in YRI population
#' @param prare_CEU proportion of causal rare SNPs in CEU population
#' @param prare_CHB proportion of causal rare SNPs in CHB population
#' @param dist distribution of effect sizes, default is "gaussian".
#'             Options can be "gaussian", "exponential", "laplace".
#'             Note that LDpred simulates non-zero effect sizes from laplace distribution.
#' @param ld linkage disequilibrium.
#'           Options can be "strong" or "weak".
#'           We simulate nearby snps together when ld is strong.
#'           We simulate snps evenly spread over all variants when ld is weak.
#' @param chunk_size the size of chunk we perform matrix multiplication
#' @param cores_used the number of cores used in parallel computing
#' @return phenotypes for three populations and effect sizes for three populations
#' @import BEDMatrix
#' @export



# https://www.r-bloggers.com/passing-arguments-to-an-r-script-from-command-lines/ set command line argument


main <- function(sim_label=1,
                 N_YRI, N_CEU, N_CHB,
                 p_YRI, p_CEU, p_CHB,
                 h2_YRI, h2_CEU, h2_CHB,
                 prare_YRI, prare_CEU, prare_CHB,
                 dist="gaussian", ld="strong",
                 chunk_size, cores_used){

  # simulate genotypes
  sim_genos_res <- sim_genos(N_YRI, N_CEU, N_CHB, sim_label)
  # extract genotype matrix
  geno_path <- sim_genos_res$genos
  geno_mat <- BEDMatrix(geno_path)
  # extract allele frequency
  af_path <- sim_genos_res$af
  af_vec <- read.delim(af_path)$ALT_FREQS
  # simulate phenotypes
  sim_phenos_res <- sim_phenos(seed=sim_label, geno_mat, chunk_size, cores_used, af_vec,
                               N_YRI, N_CEU, N_CHB,
                               p_YRI, p_CEU, p_CHB,
                               h2_YRI, h2_CEU, h2_CHB,
                               prare_YRI, prare_CEU, prare_CHB,
                               dist, ld)
  gwas_df <- perform_gwas(sim_phenos_res$phenos, sim_label, geno_path)
  return(list(phenos=sim_phenos_res$phenos,
              betas=sim_phenos_res$betas))
}







