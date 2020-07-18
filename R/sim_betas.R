#' Generate effect sizes betas
#'
#' @param M the number of snps
#' @param p the proportion of causal snps out of all snps
#' @param h2 heritability
#' @param maf_vec the vector of MAF for all snps
#' @param prare proportion of rare causal snps out of causal snps
#' @param dist distribution of effect sizes, default is "gaussian".
#'             Options can be "gaussian", "exponential", "laplace".
#'             Note that LDpred simulates non-zero effect sizes from laplace distribution.
#' @param ld linkage disequilibrium.
#'           Options can be "strong" or "weak".
#'           We simulate nearby snps together when ld is strong.
#'           We simulate snps evenly spread over all variants when ld is weak.
#' @return simulated betas of length M
#' @export
sim_betas <- function(M, p, h2, maf_vec, prare, dist, ld){
  if (dist != "gaussian" & dist != "exponential" & dist != "laplace"){
    print("Input for dist must be gaussian, or exponential, or laplace.")
  }
  if (ld != "weak" & ld != "strong"){
    print("Input for ld must be either weak or strong.")
  }
  # Initialize effect sizes
  betas <- numeric(M) * 0
  # Define rare variants as MAF < 0.01
  rare_idx <- which(maf_vec < 0.01)
  common_idx <- which(maf_vec >= 0.01)
  # The number of rare causal variants
  m_rare <- M * p * prare
  # The number of common causal variants
  m_common <- M * p - m_rare
  # Determine non-zero effect sizes based on specified ld option
  if (ld == "weak"){
    causal_idx_rare <- rare_idx[seq(m_rare) * floor(length(rare_idx)/m_rare)]
    causal_idx_common <- common_idx[seq(m_common) * floor(length(common_idx)/m_common)]
  } else if (ld == "strong"){
    causal_idx_rare <- rare_idx[seq(m_rare)]
    causal_idx_common <- common_idx[seq(m_common)]
  }
  # Sample effect sizes from the specified distribution
  if (dist == "gaussian"){
    causal_effects_rare <- rnorm(m_rare, mean=0, sd = sqrt(h2/m_rare))
    causal_effects_common <- rnorm(m_common, mean=0, sd = sqrt(h2/m_common))
  } else if (dist == "exponential"){
    causal_effects_rare <- rexp(m_rare, rate=5)
    causal_effects_common <- rnorm(m_common, rate=5)
  } else if (dist == "laplace"){
    causal_effects_rare <- rlaplace(m_rare, m=0, s = sqrt(h2/m_rare))
    causal_effects_common <- rnorm(m_common, m=0, s = sqrt(h2/m_common))
  }
  betas[causal_idx_rare] <- causal_effects_rare
  betas[causal_idx_common] <- causal_effects_common
  return(betas)
}
