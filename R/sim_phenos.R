#' Generate phenotypes for three populations YRI, CEU, and CHB.
#'
#' @param geno_mat a BEDMartix of
#'        dimension N (number of individuals from all three populations) by
#'        M (number of SNPs).
#' @param chunk_size the size of chunk we perform matrix multiplication
#' @param cores_used the number of cores used in parallel computing
#' @param maf_vec the vector of MAF for all SNPs
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
#' @return a list of phenotypes for YRI (pheno_YRI), CEU (pheno_CEU), and CHB (pheno_CHB).
#' @export

sim_phenos <- function(geno_mat, chunk_size, cores_used, maf_vec,
                      N_YRI, N_CEU, N_CHB,
                      p_YRI, p_CEU, p_CHB,
                      h2_YRI, h2_CEU, h2_CHB,
                      prare_YRI, prare_CEU, prare_CHB,
                      dist, ld){
  pheno_YRI <- sim_pheno_1pop("YRI", geno_mat, chunk_size, cores_used,
                              N_YRI, N_CEU, N_CHB,
                              p_YRI, h2_YRI, maf_vec, prare_YRI, dist, ld)
  pheno_CEU <- sim_pheno_1pop("CEU", geno_mat, chunk_size, cores_used,
                              N_YRI, N_CEU, N_CHB,
                              p_CEU, h2_CEU, maf_vec, prare_CEU, dist, ld)
  pheno_CHB <- sim_pheno_1pop("CHB", geno_mat, chunk_size, cores_used,
                              N_YRI, N_CEU, N_CHB,
                              p_CHB, h2_CHB, maf_vec, prare_CHB, dist, ld)
  return(list(
    pheno_YRI = pheno_YRI,
    pheno_CEU = pheno_CEU,
    pheno_CHB = pheno_CHB
  ))
}

