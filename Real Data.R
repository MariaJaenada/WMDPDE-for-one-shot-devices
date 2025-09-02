##################### REAL DATA EXAMPLE ########################

Temp = rep(c(308, 318, 328), each = 3)
x =matrix(1/ Temp)
tau = rep(c(10, 20, 30), times = 3) 
K = rep(10, times = 9)
observed_failures = c(3, 3, 7, 1, 5, 7, 6, 7, 9)

J = 1
I = length(tau)

theta_init = c(-10, 4000,4, -1200)
x0 = 1 / 298 # normal conditions
t0 = c(10, 20, 30)

gamma =0
gammas = seq(0, 1, 0.1)
results = numeric()

for (i in 1:I) {
  for (gamma in gammas) {
    est_theta = WMDPDE_estimator(theta_init, x, tau, K,
                                 observed_failures, gamma, method = "Nelder-Mead")
    
    a = est_theta[1:(J + 1)]
    b = est_theta[(J + 2):length(est_theta)]
    xij = c(1, x[i, ])
    
    alpha_i = exp(sum(a * xij))
    beta_i = exp(sum(b * xij))
    
    mean_life = (alpha_i * pi /beta_i) / sin(pi/ beta_i)
    
    results = rbind(results, data.frame(gamma = gamma, 
                                        theta1 = est_theta[1], 
                                        theta2 = est_theta[2], 
                                        theta3 = est_theta[3], 
                                        theta4 = est_theta[4], 
                                        alpha_i = alpha_i,
                                        beta_i = beta_i,
                                        mean_life = mean_life,
                                        R = 1 - F_model(theta = est_theta, x, tau, i),
                                        condition = i))
  }
}

head(results)

theta_init = c(-10, 4000, 4, -1200)
x0 = 1/298 # normal conditions
t0 = c(10, 20, 30)
i = 1
results = numeric()

alpha= 0.05
z = qnorm(1 - alpha/2)

for (gamma in gammas) {
  est_theta = WMDPDE_estimator(theta_init, x, tau, K,
                               observed_failures, gamma, method = "Nelder-Mead")
  
  a = est_theta[1:(J + 1)]
  b = est_theta[(J + 2):length(est_theta)]
  xij = c(1,x0) 
  
  alpha_i = exp(sum(a * xij))
  beta_i = exp(sum(b * xij))
  
  mean_life = (alpha_i * pi / beta_i) / sin(pi / beta_i)
  
  Sigma = sigma_gamma(est_theta, x, tau, I, K, sum(K), gamma)
  Sigma = Sigma /sum(K)  
  
  # standard errors and confidence intervals:
  se_theta = sqrt(diag(Sigma))
  CI_lower = est_theta - z * se_theta
  CI_upper = est_theta + z * se_theta
  
  results = rbind(results, data.frame(gamma = gamma, 
                                      theta1 = est_theta[1], 
                                      theta2 = est_theta[2], 
                                      theta3 = est_theta[3], 
                                      theta4 = est_theta[4], 
                                      alpha_i = alpha_i,
                                      beta_i = beta_i,
                                      mean_life = mean_life,
                                      R = 1 - F_model(theta = est_theta, x, tau, i),
                                      condition = i,
                                      
                                      theta1 = est_theta[1], theta1_se  = se_theta[1],
                                      theta1_lower = CI_lower[1], theta1_upper = CI_upper[1],
                                      
                                      theta2 = est_theta[2], theta2_se  = se_theta[2],
                                      theta2_lower = CI_lower[2], theta2_upper = CI_upper[2],
                                      
                                      theta3 = est_theta[3], theta3_se  = se_theta[3],
                                      theta3_lower = CI_lower[3], theta3_upper = CI_upper[3],
                                      
                                      theta4 = est_theta[4], theta4_se  = se_theta[4],
                                      theta4_lower = CI_lower[4], theta4_upper = CI_upper[4]
  )   )
  
  
}

data_condition = subset(results,condition == i)

g = ggplot(data_condition, aes(x = gamma, y =mean_life)) +
  geom_line() +
  geom_point() +
  labs(title = "Normal temperature conditions", 
       x = 'Gamma', 
       y = 'Mean lifetime') +
  theme_minimal()

print(g) 



results_gamma_0 = subset(results, gamma == 0)
print(results_gamma_0)

theta_0 = c(results_gamma_0$theta1[1],	
            results_gamma_0$theta2[1],
            results_gamma_0$theta3[1],
            results_gamma_0$theta4[1])

theoretical_prob = numeric(length(tau))

set.seed(123)

for (i in 1:length(tau)) {
  theoretical_prob[i] = F_model(theta_0,x, tau, i)
}

theoretical_prob = c(theoretical_prob, 1 - theoretical_prob) / I

sum(theoretical_prob) # equals 1

chi_squared_statistic = function(n, th) {
  N = sum(n)
  return( sum( (n - N * th)^2/(N*th) ) )
}

# Specific example MLE

print(paste("The log-logistic hypothesis should be rejected", 
            chi_squared_statistic(
              c(observed_failures, K - observed_failures), theoretical_prob) > 
              qchisq(0.95, df = 18 - (4 -1))
) 
)

print(paste("p-value of the test", 1 - pchisq(
  chi_squared_statistic(
    c(observed_failures, K - observed_failures), theoretical_prob), 
  df = 18 - (4 - 1))
) 
)

gammas = seq(0, 1, 0.1)
results_test = data.frame()

for (gamma in gammas) {
  
  est_theta = WMDPDE_estimator(theta_init, x,tau, K, 
                               observed_failures, gamma, method = "Nelder-Mead")
  
  theoretical_prob = numeric(length(tau))
  for (i in 1:length(tau)) {
    theoretical_prob[i] = F_model(est_theta, x, tau, i)  
    # We calculate p_i for each condition
  }
  
  complete_theoretical_prob = c(theoretical_prob, 1 - theoretical_prob) / I
  
  # Vector of observed frequencies (failures and successes)
  observed_freq = c(observed_failures, K - observed_failures)
  
  stat = chi_squared_statistic(observed_freq, complete_theoretical_prob)
  p_value = 1 - pchisq(stat, df = length(observed_freq) - 1)
  
  results_test = rbind(results_test, data.frame(
    gamma = gamma,
    chi_squared = stat,
    p_value = p_value
  ))
}

print(results_test)






# RESIDUAL PLOTS

gammas = c(0, 0.4)
results_test = data.frame()
results_prop = data.frame()

for (gamma in gammas) {
  
  est_theta = WMDPDE_estimator(theta_init, x, tau, K, observed_failures, gamma, method = "Nelder-Mead")
  
  theoretical_prob = numeric(length(tau))
  for (i in 1:length(tau)) {
    theoretical_prob[i] = F_model(est_theta, x, tau, i)
  }
  
  complete_theoretical_prob = c(theoretical_prob, 1 - theoretical_prob) / length(tau)
  observed_freq = c(observed_failures, K - observed_failures)
  statistic = chi_squared_statistic(observed_freq, complete_theoretical_prob)
  p_value = 1 - pchisq(statistic, df = length(observed_freq) - 1)
  
  results_test = rbind(results_test,
                       data.frame(gamma = gamma,
                                  chi_squared = statistic,
                                  p_value = p_value))
  
  # Observed vs theoretical proportions (failures)
  obs_prop = observed_failures / K
  exp_prop = theoretical_prob
  diff_prop = obs_prop - exp_prop 
  
  results_prop = rbind(results_prop,
                       data.frame(gamma = gamma,
                                  cell = 1:length(tau),
                                  Observed = obs_prop,
                                  Expected = exp_prop,
                                  Difference = diff_prop))
}

print(results_prop)
