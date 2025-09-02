source("WMDPDE")
#######################  INFLUENCE FUNCTION ####################### 


IF_gamma = function(gamma, theta0, i0, x, tau, K) {
  I = length(tau)             
  Ktotal = sum(K)             
  xij = c(1, x[i0, ])
  J =  (length(theta0) / 2) - 1
  a =  theta0[1:(J + 1)]
  b =  theta0[(J + 2):length(theta0)]
  alpha_i = exp(sum(a * xij))
  beta_i  = exp(sum(b * xij))
  F_theta = ploglogis(x = tau[i0], a = alpha_i, b = beta_i)
  R_theta = 1 - F_theta
  
  v = c(-beta_i*xij, log(tau[i0] / alpha_i)*beta_i* xij)
  
  Delta =1
  
  scalar = (K[i0] / Ktotal) *
    F_theta * R_theta *
    (F_theta^(gamma - 1) + R_theta^(gamma - 1)) *
    (F_theta - Delta)
  
  # inverse of J_gamma 
  J_inv = solve(J_gamma(theta0, I, K, x, tau, gamma))
  IF_vec = J_inv %*% v * scalar
  
  return(abs(IF_vec[1]))
}


# Influence Function vs gamma

theta0 = c(1, -0.5, 0.8, 0.4)                
I = 9
x = matrix(rep(c(0, 0.5, 1), 3), nrow = 9, ncol = 1, byrow = TRUE)
tau = c(rep(1, 3), rep(1.5, 3), rep(2.5, 3))
K = rep(100, I)

gammas = seq(0, 1, by = 0.05)

# We evaluate the IF in group i0 = 1 (the contaminated observation is in the first cell)
values_IF = sapply(gammas, IF_gamma, theta0 = theta0, i0 = 1, x = x, tau = tau, K = K)

df = data.frame(gamma = gammas, IF = values_IF)

ggplot(df, aes(x = gamma, y = IF)) + 
  geom_line() + 
  geom_point() + 
  labs(x = expression(gamma), y = "|IF|") +  
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white", color = "white"), 
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"))