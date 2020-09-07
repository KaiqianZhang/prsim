#' Helper functions for simulating phnotypes
#' Perform multiplication of a matrix by a vector in chunks
mul_chunk <- function(idx, chunk_size=100, beta, geno_mat, from, to){
  start <- from + (idx-1) * chunk_size
  end <- start + chunk_size -1
  if (end > to) end <- to
  sparse_chunk <- as(geno_mat[start:end,], "dgCMatrix")
  return(sparse_chunk %*% beta)
}

#' Helper functions for simulating phenotypes
#' Convert each dgeMatrix in the list into a vector
dge_to_vec <- function(idx, y){
  return(as.vector(y[[idx]]))
}

#' Helper functions for simulating phnotypes
#' Use multiple cores to parallely compute large matrix multiplication
mul_parallel <- function(chunk_size=100, beta, geno_mat, from, to, cores_used){
  N <- (to - from) + 1
  chunks_num <- ceiling(N/chunk_size)
  res <- mclapply(1:chunks_num, mul_chunk, chunk_size=chunk_size, beta=beta, geno_mat=geno_mat,
                  from=from, to=to, mc.cores=cores_used)
  # Each chunk of y in the whole list is a "dgeMatrix", which is hard to unlist.
  # We have to convert each chunk result back to a vector.
  res <- lapply(1:chunks_num, dge_to_vec, y=res)
  return(unlist(res))
}

#' Helper functions for simulating phnotypes
#' Simulate environmental noise
sim_E <- function(N, h2){
  epsilons <- rnorm(N, mean=0, sd=sqrt(1-h2))
  Z_epsilon <- (epsilons - mean(epsilons)) / sd(epsilons)
  E <- sqrt(1-h2) * Z_epsilon
  return(E)
}

#' Helper functions for simulating phnotypes
#' Simulate phenotypes for one specified population
sim_pheno_1pop<- function(seed, pop, geno_mat, chunk_size, cores_used,
                          N_YRI, N_CEU, N_CHB,
                          p, h2, maf_vec, prare, dist, ld){
  if (pop == "YRI"){
    N <- N_YRI
    from <- 1
    to <- N_YRI
  } else if (pop == "CEU"){
    N <- N_CEU
    from <- 1+N_YRI
    to <- N_YRI + N_CEU
  } else if (pop == "CHB"){
    N <- N_CHB
    from <- 1 + N_YRI + N_CEU
    to <- N_YRI + N_CEU + N_CHB
  } else
  {stop("The input population should be either YRI, CEU, or CHB.")}
  betas <- sim_betas(seed, dim(geno_mat)[2], p, h2, maf_vec, prare, dist, ld)
  X <- mul_parallel(chunk_size, betas, geno_mat, from, to, cores_used)
  Z_X <- (X-mean(X))/sd(X)
  G <- sqrt(h2)*Z_X
  E <- sim_E(N, h2)
  return(list(phenos = G+E, betas=betas))
}
