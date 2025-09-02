

#######################  FUNCTIONS TO GENERATE CONTAMINATION ####################### 

# For one parameter
cont_list = function(param_val, n, min_range, max_range) {
  values = seq(param_val + min_range * abs(param_val),
               param_val + max_range * abs(param_val), 
               length.out = n)
  contamination = abs((values -param_val) / param_val)
  data.frame(Contaminated_Value = values, Contamination_Level = contamination)
}

# For all parameters
generate_cont_list = function(theta, n = 10, min_range, max_range) {
  cont_list_all = lapply(theta, cont_list, n = n, min_range = min_range, max_range = max_range)
  cont_df = do.call(rbind, lapply(seq_along(cont_list_all), 
    function(i) {
    df = cont_list_all[[i]]
    df$parameter = i
    return(df)
    }
  ))
  return(cont_df)
}

# Function to generate contaminated samples
generate_contaminated_sample = function(theta_0, list, x, tau, K, seed, contaminated_cell = 1) {
  contaminated_samples = list()
  I = length(tau)
  unique_params = unique(list$parameter)
  
  for (i in unique_params) {
    sub_list = list[list$parameter == i, ]
    for (j in seq_len(nrow(sub_list))) {
      set.seed(seed)
      theta_cont = theta_0
      theta_cont[i] = sub_list$Contaminated_Value[j]
      n_cont = numeric()
      
      for (t in 1:I) {
        n_cont[t] = rbinom(1, K[t], F_model(theta_0, x, tau, t))
        if (t == contaminated_cell) {
          n_cont[t] = rbinom(1, K[t], F_model(theta_cont, x, tau, t))
        }
      }
      
      ind_sample = paste0("param_", i, "_cont_", j)
      contaminated_samples[[ind_sample]] = list(
        n_contaminated = n_cont,
        contamination_level = sub_list$Contamination_Level[j])
    }
  }
  return(contaminated_samples)
}


########################  ERROR FUNCTIONS ####################### 

compute_mae = function(real_value, estimated_value) {
  mae = mean(abs(real_value - estimated_value))
  return(mae)
}

compute_mse = function(real_value, estimated_value) {
  mse = mean((real_value - estimated_value)^2)
  return(mse)
}

compute_rmse =function(real_value, estimated_value) {
  rmse = sqrt(mean((real_value - estimated_value)^2))
  return(rmse)
}

compute_mape = function(real_value, estimated_value) {
  mape = mean(abs((real_value - estimated_value) / real_value)) * 100
  return(mape)
}


########################## SIMULATIONS ########################

# We contaminate the third parameter in the first cell
I=9  # Number of test conditions
J=1  # Number of stress factors

# Number of devices for each test condition (K_i)
K = rep(100,I) 

# Initial parameters 
theta_init = c(1, 0, 1, 0)

# Stress factors (x_ij)
x = matrix(rep(c(0, 0.5, 1), 3), nrow = 9, ncol= 1, byrow = TRUE)
tau = c(rep(1, 3), rep(1.5, 3), rep(2.5,3))
theta_0 = c(1, -0.5, 0.8, 0.4) # true value to estimate

results = data.frame(
  parameter = numeric(),
  Error = numeric(),
  Contamination_Level = numeric(),
  Gamma = numeric(),
  Theta_1 = numeric(),
  Theta_2 = numeric(),
  Theta_3 = numeric(),
  Theta_4 = numeric(),
  Level = numeric()
)

CONTAMINATED_CELL =1

# Wald test
m = matrix(c(0, 0, 1, 0), 4, 1)
M = matrix(0, 4, 1)
M[3, 1] = 1

theta_0 = c(1, -0.5, 0.8, 0.4)
theta_0_H1 = c(1, -0.5, 0.35, 0.4)

n_sim = 3000
n = 10
alpha = 0.05
min_range =-0.5
max_range = 0

contamination_list = generate_cont_list(theta_0, n, min_range, max_range)

gamma_values = seq(0, 1, by = 0.2)

results = data.frame(
  parameter = numeric(),
  Error = numeric(),
  Contamination_Level = numeric(),
  Gamma = numeric(),
  Theta_1 = numeric(),
  Theta_2 = numeric(),
  Theta_3 = numeric(),
  Theta_4 = numeric(),
  Level = numeric()
)

contamination_list =generate_cont_list(theta_0, n, min_range, max_range)

CONTAMINATED_CELL = 1

for (gamma in gamma_values) {
  sim_errors = data.frame(Cont_Level = numeric(), 
                          Error_1 = numeric(), Error_2 = numeric(), 
                          Error_3 = numeric(), Error_4 = numeric(), 
                          Theta_1 = numeric(), Theta_2 = numeric(), 
                          Theta_3 = numeric(), Theta_4 = numeric(),
                          n_contaminated_1 = numeric(), n_contaminated_2 = numeric(), 
                          n_contaminated_3 = numeric(), n_contaminated_4 = numeric(), 
                          n_contaminated_5 = numeric(), n_contaminated_6 = numeric(), 
                          n_contaminated_7 = numeric(), n_contaminated_8 = numeric(), 
                          n_contaminated_9 = numeric())
  
  levels_df = data.frame()
  
  for (sim in 1:n_sim) {
    # we ensure that for each gamma the samples are the same:
    seed = 1234 + sim
    sub_list = contamination_list[contamination_list$parameter == 3, ]  
    contaminated_samples = generate_contaminated_sample(theta_0, sub_list, x, tau, K, seed, contaminated_cell = CONTAMINATED_CELL)
    
    for (sample_idx in seq_along(contaminated_samples)) {
      theta_contaminated = theta_0
      theta_contaminated[3] = sub_list$Contaminated_Value[sample_idx]  
      
      sample =contaminated_samples[[sample_idx]]
      
      estimates = WMDPDE_estimator(theta_init, x,
                                   tau, K, sample$n_contaminated, 
                                   gamma, method = "Nelder-Mead")
      
      error_1 = compute_rmse(theta_0[1], estimates[1])
      error_2 = compute_rmse(theta_0[2], estimates[2])
      error_3 = compute_rmse(theta_0[3], estimates[3])
      error_4 = compute_rmse(theta_0[4], estimates[4])
      error = compute_rmse(theta_0, estimates)
      
      sim_errors = rbind(sim_errors, data.frame(
        Cont_Level = sub_list$Contamination_Level[sample_idx],
        Error_1 = error_1,Error_2 = error_2, 
        Error_3 = error_3, Error_4 = error_4, Error = error,
        Theta_1 = estimates[1], Theta_2 = estimates[2],
        Theta_3 = estimates[3], Theta_4 = estimates[4],
        n_contaminated_1 = sample$n_contaminated[1], 
        n_contaminated_2 = sample$n_contaminated[2],
        n_contaminated_3 = sample$n_contaminated[3], 
        n_contaminated_4 = sample$n_contaminated[4],
        n_contaminated_5 = sample$n_contaminated[5], 
        n_contaminated_6 = sample$n_contaminated[6],
        n_contaminated_7 = sample$n_contaminated[7], 
        n_contaminated_8 = sample$n_contaminated[8],
        n_contaminated_9 = sample$n_contaminated[9]
      ))
      
      sigma = sigma_gamma(estimates, x, tau, I, K, sample$n_contaminated, gamma)
      
      if (is.null(sigma)) {
        sim = sim + 1
        count = count + 1
      } else {
        t0 = sum(K) * (estimates - theta_0) %*% m %*% solve(t(M) %*% sigma %*% M) %*% t(t((estimates - theta_0)) %*% m)
        t0_H1 = sum(K) * (estimates - theta_0_H1) %*% m %*% solve(t(M) %*% sigma %*% M) %*% t(t((estimates - theta_0_H1)) %*% m)
        
        level_val = as.numeric(t0 > qchisq(1 - alpha, df = 1)) 
        # convert TRUE/FALSE to 1/0
        power_val = as.numeric(t0_H1 > qchisq(1 - alpha, df = 1))
        
        levels_df = rbind(levels_df, data.frame(
          level =level_val,
          Cont_Level = sub_list$Contamination_Level[sample_idx],
          power = power_val
        ))
      }
    }
  }
  
  avg_errors = sim_errors %>%
    group_by(Cont_Level) %>%
    summarise(
      Mean_Error_1 = mean(Error_1), Mean_Error_2 = mean(Error_2),
      Mean_Error_3 = mean(Error_3), Mean_Error_4 = mean(Error_4),
      Mean_Error = mean(Error),
      Theta_1 = mean(Theta_1), Theta_2 = mean(Theta_2),
      Theta_3 = mean(Theta_3), Theta_4 = mean(Theta_4),
      n_contaminated_1 = mean(n_contaminated_1), n_contaminated_2 = mean(n_contaminated_2),
      n_contaminated_3 = mean(n_contaminated_3),n_contaminated_4 = mean(n_contaminated_4),
      n_contaminated_5 = mean(n_contaminated_5), n_contaminated_6 = mean(n_contaminated_6),
      n_contaminated_7 = mean(n_contaminated_7), n_contaminated_8 = mean(n_contaminated_8),
      n_contaminated_9 = mean(n_contaminated_9)
    )
  
  # we compute test significance level and power
  
  avg_levels = levels_df %>%
    group_by(Cont_Level) %>%
    summarise(Mean_Level = mean(level), Mean_Power = mean(power))
  
  results = rbind(results, 
                  cbind(
                    data.frame(
                      Gamma = gamma,
                      Error_1 = avg_errors$Mean_Error_1,
                      Error_2= avg_errors$Mean_Error_2,
                      Error_3 = avg_errors$Mean_Error_3,
                      Error_4 = avg_errors$Mean_Error_4,
                      Error = avg_errors$Mean_Error,
                      Cont_Level = avg_errors$Cont_Level,
                      Level = avg_levels$Mean_Level,
                      Power = avg_levels$Mean_Power),
                    avg_errors[, c("Theta_1", "Theta_2", 
                                   "Theta_3", "Theta_4",
                                   "n_contaminated_1", "n_contaminated_2", 
                                   "n_contaminated_3", "n_contaminated_4",
                                   "n_contaminated_5", "n_contaminated_6", 
                                   "n_contaminated_7", "n_contaminated_8",
                                   "n_contaminated_9")]
                  )
  )
}

