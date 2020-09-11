# test file
library(prsim)

res <- main(sim_label=2,
            N_YRI=200, N_CEU=0, N_CHB=0,
            p_YRI=0.3, p_CEU=0, p_CHB=0,
            h2_YRI=0.5, h2_CEU=0, h2_CHB=0,
            prare_YRI=0.1, prare_CEU=0, prare_CHB=0,
            dist="gaussian", ld="strong",
            chunk_size=100, cores_used=3)
