
## Set up immigration-death model

# Data
id_ods = readRDS(file = "data/id.rds")
# Prior
id_log_prior = function(param)  sum(dnorm(param, mean = 0, sd = 10, log = TRUE))
# Initial parameter vector
pars_init = c(4.0, 0.8)

# SDE 
id_sde = function(theta) {
  theta = theta
  drift = function(x) theta[1] - theta[2] * x
  diffusion = function(x) sqrt(theta[1] + theta[2] * x)
  list(drift = drift, diffusion = diffusion)
}

plot_theta = function(theta) {
  op = par(mar = c(3,3,2,1),
          mgp = c(2,0.4,0), tck=-.01,
          cex.axis = 0.9, las=1, mfrow=c(ncol(theta),3))
  on.exit(par(op))
  i = 1
  for(i in 1:ncol(theta)) {
    par_lab = paste0("log(c", i, ")")
    plot(theta[,i], ylab = par_lab, xlab = "Iteration", type = "l")
    lines(lowess(theta[,i]), col= "steelblue", lwd = 1.5)
    acf(theta[,i], main = NA)
    plot(density(theta[,i]), xlab=par_lab, main = NA, 
         ylab = "Density")
  }
}
