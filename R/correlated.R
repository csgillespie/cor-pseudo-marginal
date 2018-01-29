source("R/helper.R")

#As bridgeSimQ but supply Brownian innovations
bridgeInc = function(end_points, deltat, inter, sim, bmvec) {
  #Supply bm increments, elements Wti-Wti-1 are N(0,deltat) distributed
  #bmvec is a (n-1)-vector of N(0,deltat) quantities
  n = inter/deltat
  vec = numeric(n+1)
  vec[c(1, n + 1)] = end_points
  for(i in 1:(n-1)) {
    tm = inter - (i-1) * deltat
    t = inter - i * deltat
    vec[i + 1] = (vec[i] + (end_points[2] - vec[i])* deltat/tm) +
      sim$diffusion(vec[i]) * bmvec[i] * sqrt(t/tm)
    vec[i + 1] = max(vec[i + 1], 1e-5)
  }
  vec
}

#Importance sampler to estimate transition density p(xT|x0,theta)
imp_samp = function(N, end_points, deltat, inter, sim, bmmat) {
  #bmmat is an N * (n-1) matrix of N(0,deltat) quantities
  wts = numeric(N)
  for(i in 1:N) {
    mat = bridgeInc(end_points, deltat, inter, sim, bmmat[i,,])
    wts[i] = exp(EMdensity(mat, deltat, inter, sim) -
                   Qdensity(mat,deltat, inter, sim))
  }
  mean(wts)
}


         
correlated = function(sde, obs, sigma, log_prior, pars_init,
                   deltat, inter, 
                   iters = 1000, N = 1, rho = 0) {
  n=length(obs) # no. obs
  n2 = (inter/deltat)-1 #no. BM increments in each interval
  bmarray_cur = array(rnorm(N*n2*(n-1), 0, sqrt(deltat)),
                      dim = c(N, n2, n-1))
  p = length(pars_init)  #no. params
  mat = matrix(0, ncol = p + 1, nrow = iters)
  cur = log(pars_init) #current param value (log-scale)
  ll_cur = -Inf #log-likelihood of current value (to accept first candidate)
  mat[1,] = c(cur, 0) #params on log-scale and candidate log-like
  count = 0

  for (i in 2:iters) {
    #update bm increments
    bmarray_prop = rho * bmarray_cur +
      sqrt(1-rho^2) * rnorm(N * n2 * (n-1), 0, sqrt(deltat))

    prop = rmvn(cur, sigma)
    sim_prop = sde(exp(prop))
    ll_prop = 0
    for(j in 1:(n-1)) {
      like = imp_samp(N, end_points = obs[j:(j+1)],
                      deltat, inter, sim_prop, bmarray_prop[,,j, drop = FALSE])
      ll_prop = ll_prop + log(like)
    }

    #Evaluate (log) acceptance probability
    accept_prob = ll_prop + log_prior(prop) -
      ll_cur - log_prior(cur)

    if (log(runif(1)) < accept_prob) {
      cur = prop
      ll_cur = ll_prop
      bmarray_cur = bmarray_prop
      count = count+1  #track no. acceptances
    }
    mat[i,] = c(cur, ll_cur)
  }
  message(signif(count/(iters-1), 3)*100)
  return(mat[2:iters, -ncol(mat)]) #return (log) parameter samples
}
