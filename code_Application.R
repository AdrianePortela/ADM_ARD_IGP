# =========================================================
# Application
# =========================================================

load("wind_turbine_dat.Rdata")

library(dplyr)

estimacao <- function(wind_turbine_data) {
  
  # 1. Identify maintenance positions (adapted for the specific columns)
  identify_repeated_positions <- function(data, unit_id) {
    unit_data <- data %>% filter(Objeto == unit_id)
    time_vector <- unit_data$Time
    
    frequency <- table(time_vector)
    repeated_values <- names(frequency[frequency > 1])
    
    repeated_positions <- unlist(lapply(repeated_values, function(value) {
      which(time_vector == value)
    }))
    
    repeated_positions <- sort(repeated_positions)
    positions <- repeated_positions[seq_along(repeated_positions) %% 2 != 0]
    
    return(positions)
  }
  
  # 2. Compute the required quantities for the likelihood estimation
  # Loop iterates through the levels of the 's' factor (acceleration level)
  s_levels <- unique(wind_turbine_data$s)
  
  all_Z_increments_global <- c()
  all_time_increments_global <- c()
  s_global <- NULL
  
  for (j in 1:length(s_levels)) {
    s_value <- s_levels[j]
    
    # Filter data for the current acceleration level
    data_filtered <- wind_turbine_data %>% filter(s == s_value)
    
    # Identify positions using the first object of this filtered block
    primeiro_objeto <- data_filtered$Objeto[1]
    positions <- identify_repeated_positions(data_filtered, primeiro_objeto)
    
    # Number of units at this acceleration level
    m <- length(table(data_filtered$Objeto))
    objetos_no_nivel <- unique(data_filtered$Objeto)
    
    for (i in 1:m) {
      subset_data <- subset(data_filtered, Objeto == objetos_no_nivel[i])
      
      Z <- subset_data$Y
      t <- subset_data$Time
      
      # Compute increments
      Z_increment <- Z[-1] - Z[-length(Z)]
      time_increment <- t[-1] - t[-length(t)]
      
      # Remove maintenance positions
      Z_increment_clean <- Z_increment[-c(positions)]
      time_increment_clean <- time_increment[-c(positions)]
      
      # Update global vectors
      all_Z_increments_global <- c(all_Z_increments_global, Z_increment_clean)
      all_time_increments_global <- c(all_time_increments_global, time_increment_clean)
    }
    
    s_rep <- rep(s_value, length(all_Z_increments_global) - length(s_global))
    s_global <- c(s_global, s_rep)
  }
  
  # Structure the final data frame for estimation
  estimation_data <- data.frame(
    increment = all_Z_increments_global,
    dt = all_time_increments_global,
    A = s_global
  )
  
  # 3. Negative Log-Likelihood function
  LogLik <- function(par, data) {
    delta.t <- par[1] * data$dt * exp(par[3] * data$A)
    
    mu <- par[1]
    lambda <- par[2]
    beta <- par[3]
    
    increment <- data$increment
    
    t1 <- log(sqrt(lambda) * delta.t)
    t2 <- -(lambda * (increment - delta.t)^2) / (2 * increment)
    
    l <- -(sum(t1) + sum(t2))
    
    if (!is.finite(l)) l <- 1e50
    
    return(l %>% as.numeric())
  }
  
  # 4. Optimization process with tryCatch
  start_values <- c(1, 1, 0)
  method <- "BFGS"
  
  mod <- tryCatch(
    optim(
      fn = LogLik,
      par = start_values,
      hessian = TRUE,
      method = method,
      data = estimation_data
    ),
    error = function(e) {e}
  )
  
  # 5. Extract results and compute output metrics
  e.mu <- mod$par[1]
  e.lambda <- mod$par[2]
  e.beta <- mod$par[3]
  
  mcov <- solve(mod$hessian)
  
  sd.mu <- sqrt(mcov[1,1])
  sd.lambda <- sqrt(mcov[2,2])
  sd.beta <- sqrt(mcov[3,3])
  
  IC_mu <- e.mu + c(-1,1) * qnorm(0.975) * sd.mu
  IC_lambda <- e.lambda + c(-1,1) * qnorm(0.975) * sd.lambda
  IC_beta <- e.beta + c(-1,1) * qnorm(0.3) * sd.beta # Keeping qnorm(0.975) as standard for 95% CI
  IC_beta <- e.beta + c(-1,1) * qnorm(0.975) * sd.beta 
  
  
  # Log-likelihood and information criteria
  logLik <- -mod$value
  
  n <- nrow(estimation_data)
  p <- 3
  
  AIC <- 2 * p - 2 * logLik
  BIC <- p * log(n) - 2 * logLik
  
  results <- c(
    e.mu, sd.mu, IC_mu[1], IC_mu[2], IC_mu[2] - IC_mu[1],
    e.lambda, sd.lambda, IC_lambda[1], IC_lambda[2], IC_lambda[2] - IC_lambda[1],
    e.beta, sd.beta, IC_beta[1], IC_beta[2], IC_beta[2] - IC_beta[1],
    logLik, AIC, BIC
  )
  
  results <- results %>%
    setNames(
      c(
        "emv_mu", "sd_mu", "CI_lower_mu", "CI_upper_mu", "CI_width_mu",
        "emv_lambda", "sd_lambda", "CI_lower_lambda", "CI_upper_lambda", "CI_width_lambda",
        "emv_beta", "sd_beta", "CI_lower_beta", "CI_upper_beta", "CI_width_beta",
        "logLik", "AIC", "BIC"
      )
    )
  
  return(results)
}

final_results <- estimacao(wind_turbine_data)
print(resultados)
