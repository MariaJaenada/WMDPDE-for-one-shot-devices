
library(VaRES)
library(ggplot2)
library(reshape2)
library(dplyr)
library(parallel)
library(nleqslv)
library(gridExtra)
library(optimx)


######################## DISTRIBUTION FUNCTION ####################### 

F_model = function(theta, x, tau, i) {
  # We extract parameters a and b from vector theta = (a,b)
  J = floor(length(theta) / 2) - 1
  a = theta[1:(J + 1)]
  b = theta[(J + 2):length(theta)]
  
  # We add a 1 at the beginning of vector x for the intercept
  xij = c(1, x[i, ])
  
  alpha_i = exp(sum(a * xij))
  beta_i= exp(sum(b * xij))
  F_theta = ploglogis(x = tau[i], a = alpha_i, b = beta_i)
  
  return(F_theta)
}

########################  DISTANCE ####################### 

d = function(theta, x, tau, K, n, gamma) {
  
  K_total = sum(K)
  I =length(tau)
  F_theta = sapply(1:I, function(i) F_model(theta, x, tau, i))
  
  if (gamma == 0) {
    term1 = -(n / K) * log(F_theta)
    term2 = -((K - n) / K) * log(1 - F_theta)
    dist = sum(term1 + term2)
  } else {
    term1 = F_theta^(1 + gamma) + (1 - F_theta)^(1 + gamma)
    term2 = (F_theta^gamma * n/K) + ((1 - F_theta)^gamma * (K - n)/K)
    dist =sum(K * (term1 - (1 + 1/gamma) * term2))
  }
  
  return(dist/K_total)
}

########################  FUNCTION TO COMPUTE THE ESTIMATOR ####################### 

WMDPDE_estimator = function(theta_init, x, tau, K, n, gamma, method) {
  result = optim(
    par = theta_init,
    f = d,
    x = x,
    tau = tau,
    K = K,
    n = n,
    gamma = gamma,
    method = method,
    control = list(fnscale = 1, reltol = 1.0e-14) 
  )
  return(result$par)
}



########################  MATRICES TO CALCULATE SIGMA ####################### 

M_i = function(x_i, tau_i, alpha_i, beta_i) {
  x_i = as.matrix(x_i)
  term1 =log(tau_i / alpha_i)
  term2 = beta_i^2
  M = rbind(cbind(term2 * (x_i %*% t(x_i)), 
                  -term1 * term2 * (x_i %*% t(x_i))),
            cbind(-term1 * term2 * (x_i %*% t(x_i)), 
                  (term1 * beta_i)^2 * (x_i %*% t(x_i))))
  return(M)
}

# J_gamma(theta)
J_gamma = function(theta, I, K, x, tau, gamma) {
  J = (length(theta) / 2) - 1
  J_matrix = 0
  K_total = sum(K)
  
  for (i in 1:I) {
    R_theta = 1 - F_model(theta, x, tau, i)
    a = theta[1:(J + 1)]
    b =theta[(J + 2):length(theta)]
    xij = c(1, x[i, ])
    alpha_i = exp(sum(a * xij))
    beta_i = exp(sum(b * xij))
    M = M_i(xij, tau[i], alpha_i, beta_i)
    term = (R_theta * F_model(theta, x, tau, i))^2 *
      (F_model(theta, x, tau, i)^(gamma - 1) + R_theta^(gamma - 1))
    J_matrix = J_matrix + (K[i] / K_total) * M * term
  }
  return(J_matrix)
}

# K_gamma(theta)
K_gamma = function(theta, I, K, x, tau, gamma) {
  J = (length(theta) / 2) - 1
  K_matrix = 0
  K_total = sum(K)
  
  for (i in 1:I) {
    R_theta = 1 - F_model(theta, x, tau, i)
    a = theta[1:(J + 1)]
    b = theta[(J + 2):length(theta)]
    xij = c(1, x[i, ])
    alpha_i = exp(sum(a * xij))
    beta_i = exp(sum(b * xij))
    M = M_i(xij, tau[i], alpha_i, beta_i)
    term = (F_model(theta, x, tau, i)^(gamma - 1) + R_theta^(gamma - 1))^2 *
      F_model(theta, x, tau, i)^3 * R_theta^3
    K_matrix =K_matrix + (K[i]/ K_total) * M * term
  }
  return(K_matrix)
}

# Sigma_gamma(theta)
sigma_gamma = function(theta, x, tau, I, K, n, gamma) {
  J_inv = tryCatch(solve(J_gamma(theta, I, K, x, tau, gamma)), error = function(e) NULL)
  if (is.null(J_inv)) {
    sigma = NULL
  } else {
    sigma = J_inv %*% K_gamma(theta, I, K, x, tau, gamma) %*% J_inv
    if (is.nan(sigma[1, 1])) sigma = NULL
  }
  return(sigma)
}