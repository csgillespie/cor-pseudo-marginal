# Set-up
library("coda")
source("R/id_set_up.R")

## Modified bridge scheme
source("R/mis.R")

################################################
# Optimal tuning matrix
################################################
mis_sigma = matrix(c(0.012, 0.0052, 0.0052, 0.0046),
                ncol = 2, byrow = TRUE)

# Run sampler & time
system.time({id_mis <- mis(sde = id_sde, obs = id_ods, 
                           sigma = mis_sigma, pars_init = pars_init,
                           log_prior = id_log_prior, deltat = 0.2, inter = 1, 
                           iters = 20000)})

# Plot output
plot_theta(id_mis)
effectiveSize(id_mis)

################################################
#Correlated pseudo-marginal Metropolis-Hastings#
################################################

#Better tuning matrix:
source("R/correlated.R")
cor_sigma = matrix(c(0.023, 0.0061, 0.0061, 0.0053),
                ncol = 2, byrow = TRUE)

system.time(id_cor <- correlated(sde = id_sde, obs = id_ods,
                                 sigma = cor_sigma, pars_init = pars_init,
                                 log_prior = id_log_prior, deltat = 0.2, inter = 1, 
                                 iters = 20000, N = 1, rho = 0.99))

plot_theta(id_cor)
effectiveSize(id_cor)
