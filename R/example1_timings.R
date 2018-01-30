# Set-up
library("coda")
source("R/id_set_up.R")
reps = 2 # replications; average over timings

## Modified bridge scheme
source("R/mis.R")

################################################
# Optimal tuning matrix
################################################
mis_sigma = matrix(c(0.012, 0.0052, 0.0052, 0.0046),
                   ncol = 2, byrow = TRUE)

# Run sampler & time
dd_mis = data.frame(elapsed = 0, ess = numeric(reps), ess_sec = 0)
for(i in 1:reps) {
  tim = system.time({id_mis <- mis(sde = id_sde, obs = id_ods, 
                                   sigma = mis_sigma, pars_init = pars_init,
                                   log_prior = id_log_prior, deltat = 0.2, inter = 1, 
                                   iters = 20000)})
  dd_mis[i,] = c(tim["elapsed"], min(effectiveSize(id_mis)), min(effectiveSize(id_mis))/tim["elapsed"])
}
signif(colMeans(dd_mis),3)
saveRDS(dd_mis, file = "data/mis_timings.rds")
#16.11 17.73


################################################
#Correlated pseudo-marginal Metropolis-Hastings#
################################################

#Better tuning matrix:
source("R/correlated.R")
cor_sigma = matrix(c(0.023, 0.0061, 0.0061, 0.0053),
                   ncol = 2, byrow = TRUE)

#Run sampler and plot output appropriately
rhos = c(0, 0.5, 0.9, 0.95,0.99)
dd_cor = data.frame(rho = 0, elapsed = 0, ess = numeric(length(rhos)), N=0)
i = 1
k = 1
for(r in 1:reps) {
  for(i in seq_along(rhos)) {
    for(n in c(1, 2, 10, 100)) {
      t = system.time({id_cor <- correlated(sde = id_sde, obs = id_ods,
                                           sigma = cor_sigma, pars_init = pars_init,
                                           log_prior = id_log_prior, deltat = 0.2, inter = 1, 
                                           iters = 20000, N = n, rho = rhos[i])})
                      
                      dd_cor[k,] = c(rhos[i], t[3], ess = min(effectiveSize(id_cor)), N=n)
                      k = k + 1
    }
  }
}
saveRDS(dd_cor, file= "data/cor_timings.rds")
#dd = readRDS("res.rds")

library("dplyr")
library("tidyr")


sum = dd_cor %>%
  group_by(N, rho) %>%
  summarise(time = mean(elapsed), ess = mean(ess), ess_sec = mean(ess)/mean(elapsed))

## PMMH
dd_pmmh = data.frame(elapsed = 0, ess = numeric(reps), ess_sec = 0)
for(i in 1:reps) {
  tim = system.time({id_pmmh <- correlated(sde = id_sde, obs = id_ods,
                                           sigma = cor_sigma, pars_init = pars_init,
                                           log_prior = id_log_prior, deltat = 0.2, inter = 1, 
                                           iters = 20000, N = 50, rho = 0)})
  dd_mis[i,] = c(tim["elapsed"], min(effectiveSize(id_pmmh)), min(effectiveSize(id_pmmh))/tim["elapsed"])
}
#saveRDS(dd_mis, file = "data/pmmh_timings.rds")

signif(colMeans(dd_mis), 3)


RPushbullet::pbPost("note", title = "Completed")
  
