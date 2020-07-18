#' Generate effect sizes betas
#'
#' @param M the number of snps
#' @param m the number of causal snps
#' @param h2 heritability
#' @param maf_vec the vector of MAF for all snps
#' @param prare proportion of rare causal snps
#' @return simulated betas of length M
#' @export
sim_betas <- function(M, m, h2, maf_vec, prare){
  betas <- numeric(M) * 0
  # Define rare variants as MAF < 0.01
  rare_idx <- which(maf_vec < 0.01)
  common_idx <- which(maf_vec >= 0.01)
  # The number of rare causal variants
  m_rare <- m * prare
  causal_idx_rare <- rare_idx[seq(m_rare) * floor(length(rare_idx)/m_rare)]
  causal_effects_rare <- rnorm(m_rare, mean=0, sd = sqrt(h2/m_rare))
  betas[causal_idx_rare] <- causal_effects_rare
  # The number of common causal variants
  m_common <- m - m_rare
  causal_idx_common <- common_idx[seq(m_common) * floor(length(common_idx)/m_common)]
  causal_effects_common <- rnorm(m_common, mean=0, sd = sqrt(h2/m_common))
  betas[causal_idx_common] <- causal_effects_common
  return(betas)
}
