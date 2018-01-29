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
mis_timings = numeric(reps)
for(i in 1:reps) {
  tim = system.time({id_mis <- mis(sde = id_sde, obs = id_ods, 
                                   sigma = mis_sigma, pars_init = pars_init,
                                   log_prior = id_log_prior, deltat = 0.2, inter = 1, 
                                   iters = 20000)})
  mis_timings[i] = min(effectiveSize(id_mis))/tim["elapsed"]
}
mis_timings
#16.11 17.73


################################################
#Correlated pseudo-marginal Metropolis-Hastings#
################################################

#Better tuning matrix:
source("R/correlated.R")
cor_sigma = matrix(c(0.023, 0.0061, 0.0061, 0.0053),
                   ncol = 2, byrow = TRUE)

#Run sampler and plot output appropriately
rhos = c(0, 0.5, 0.9, 0.99)
dd = data.frame(rho = 0, elapsed = 0, ess = numeric(length(rhos)), N=0)
i = 1
k = 1
for(r in 1:2) {
  for(i in seq_along(rhos)) {
    for(n in c(1, 2, 10, 100)) {
      t = system.time({id_cor <- correlated(sde = id_sde, obs = id_ods,
                                           sigma = cor_sigma, pars_init = pars_init,
                                           log_prior = id_log_prior, deltat = 0.2, inter = 1, 
                                           iters = 20000, N = n, rho = rhos[i])})
                      
                      dd[k,] = c(rhos[i], t[3], ess = min(effectiveSize(id_cor)), N=n)
                      k = k + 1
    }
  }
}
#saveRDS(dd, file= "res.rds")
#dd = readRDS("res.rds")

library("dplyr")
library("tidyr")
head(dd)

sum = dd %>%
  group_by(N, rho) %>%
  summarise(ess_sec = mean(ess)/mean(elapsed), ess = mean(ess)) %>%
  gather(type, value, -N, -rho)
sum

