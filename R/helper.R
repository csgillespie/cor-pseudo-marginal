########################################
#Helper functions for inference schemes#
########################################
#Function to generate a single draw of a p-variate normal with mean m and covariance V
rmvn = function(m, V) {
  p = length(m)
  z = rnorm(p)
  return(m + t(chol(V)) %*% z)
}

#log of joint density of (univariate) path x under the Euler-Maruyama approximation
EMdensity = function(x, deltat, inter, sim) {
  n = inter/deltat
  x_n = x[1:n]
  means = x_n + sim$drift(x_n) * deltat
  sds = sim$diffusion(x_n) * sqrt(deltat)
  log_prob = dnorm(x[2:(n + 1)], means, sds, log = TRUE)
  sum(log_prob)
}

#log of joint density of (univariate) bridge x using modified diffusion bridge construct
Qdensity = function(x, deltat, inter, sim) {
  n = inter/deltat
  i = seq_len(n - 1)
  t = i * deltat
  tm = (i-1) * deltat

  means = x[i] + (x[n+1] - x[i]) / (inter - tm) * deltat
  sds = sim$diffusion(x[i]) * sqrt(deltat) * sqrt((inter - t) / (inter - tm))
  sum(dnorm(x[i+1], means, sds, log = TRUE))
}

#Generate a single realisation from the modified diffusion bridge construct
bridgeSimQ = function(end_points, deltat, inter, sim) {
  n = inter/deltat
  vec = numeric(n + 1)
  vec[c(1, n + 1)] = end_points

  for(i in 1:(n-1)) {
    tm = (i-1)*deltat
    t = i*deltat
    vec[i+1] = rnorm(1,
                     vec[i] + (end_points[2] - vec[i]) * deltat / (inter - tm),
                     sim$diffusion(vec[i]) * sqrt(deltat) * sqrt((inter - t) / (inter - tm)))
    #process is non-negative but EM is not
    vec[i+1] = max(vec[i+1], 1e-5)
  }
  vec
}
