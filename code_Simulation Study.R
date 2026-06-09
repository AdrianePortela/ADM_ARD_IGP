# =========================================================
# Simulation Study (SimDesign) - Full Metrics 
# Compatible with Toy Example
# =========================================================
require(statmod)
require(stats)
require(graphics)
require(tidyverse)
require(gridExtra)
library(SimDesign)
require(tibble)
library(ggplot2)
library(ggh4x)
library(scales)


# =========================================================
# Design
# =========================================================

# Simulation design: combinations of parameter values considered in the study
Design <- createDesign(
  t_max = c(4), # or outhers 
  m = c(50, 100, 200),                 # number of units
  nj = c(0, 3, 5),
  k = c(2, 3, 4),
  mu = c(2, 5, 7),
  lambda = c(10, 5, 1),
  beta = c(-1, -2, -4) # and  c(1, 2, 4)
)

#condition<-Design[1,] #block test

# =========================================================
# Generate 
# =========================================================
Generate <- function(condition, fixed_objects = NULL) {
  
  # Parameters
  t_max  <- condition$t_max
  nj     <- condition$nj
  k      <- condition$k
  mu     <- condition$mu
  lambda <- condition$lambda
  m      <- condition$m
  beta   <- condition$beta
  
  theta <- c(0.3, 0.5, 0.7, 0.9)  # maintenance effects
  A     <- c(0, 0.5, 0.75)         # acceleration levels
  
  
  # =========================================================
  # Generate one degradation trajectory
  # =========================================================
  generate_trajectory <- function(t_max, nj, k, mu, lambda,
                                  theta, id = 1, beta, A) {
    
    N <- nj * (k + 1)
    n_points <- N + k + 2
    
    time_x <- seq(0, t_max, length.out = n_points)
    dt <- t_max / (n_points - 1)
    
    # Latent process X(t)
    increments <- rinvgauss(
      n_points - 1,
      mean = mu * dt * exp(beta * A),
      shape = (mu * dt * exp(beta * A))^2 * lambda
    )
    
    X <- c(0, cumsum(increments))
    
    
    # Maintenance times
    t_change <- seq(0, t_max, t_max / (k + 1))
    t_MP <- t_change[-c(1, length(t_change))]
    
    time_z <- sort(c(time_x, t_MP))
    
    
    # Observed process Z(t)
    Z <- numeric(length(time_z))
    
    j <- length(Z) / (k + 1)
    l <- (length(X) - 1) / (k + 1)
    
    for (i in 1:k) {
      
      Z[1:j] <- X[1:(l + 1)]
      
      reduction <- 0
      
      for (p in 1:i) {
        reduction <- reduction +
          theta[p] * (X[(p * l + 1)] - X[((p - 1) * l + 1)])
      }
      
      Z[(i * j + 1):((i + 1) * j)] <-
        X[(i * l + 1):((i + 1) * l + 1)] - reduction
    }
    
    
    # Replicate X values at maintenance times (for plotting)
    X_plot <- numeric(length(Z))
    
    for (i in 1:k) {
      X_plot[1:j] <- X[1:(l + 1)]
      X_plot[(i * j + 2):((i + 1) * j)] <-
        X[(i * l + 2):((i + 1) * l + 1)]
    }
    
    for (i in 1:k) {
      X_plot[(i * j + 1)] <- X[(i * l + 1)]
    }
    
    
    # Output data frames
    df_x <- data.frame(
      id = rep(paste0("OBJ_", id)),
      time = time_x,
      X = X
    )
    
    df_z <- data.frame(
      id = rep(paste0("OBJ_", id)),
      time = time_z,
      Z = Z,
      X = X_plot
    )
    
    return(
      list(
        df_x = df_x,
        df_z = df_z,
        t_MP = t_MP
      )
    )
  }
  
  
  # =========================================================
  # Generate complete dataset
  # =========================================================
  generate_dataset <- function(t_max, nj, k, mu, lambda,
                               theta, m, beta, A) {
    
    dataset_list <- list()
    
    t_change <- seq(0, t_max, t_max / (k + 1))
    t_MP <- t_change[-c(1, length(t_change))]
    
    for (j in 1:length(A)) {
      
      list_df_x <- list()
      list_df_z <- list()
      
      for (i in 1:m) {
        
        trajectory <- generate_trajectory(
          t_max = t_max,
          nj = nj,
          k = k,
          mu = mu,
          lambda = lambda,
          theta = theta,
          id = i,
          beta = beta,
          A = A[j]
        )
        
        trajectory$df_x$A <- A[j]
        trajectory$df_z$A <- A[j]
        
        list_df_x[[i]] <- trajectory$df_x
        list_df_z[[i]] <- trajectory$df_z
      }
      
      df_x_A <- do.call(rbind, list_df_x)
      df_z_A <- do.call(rbind, list_df_z)
      
      dataset_list[[j]] <- list(
        df_x = df_x_A,
        df_z = df_z_A,
        A = A[j]
      )
    }
    
    
    df_x_final <- do.call(
      rbind,
      lapply(dataset_list, function(x) x$df_x)
    )
    
    df_z_final <- do.call(
      rbind,
      lapply(dataset_list, function(x) x$df_z)
    )
    
    
    return(
      list(
        df_x = df_x_final,
        df_z = df_z_final,
        t_MP = t_MP,
        theta = theta
      )
    )
  }
  
  
  dat <- generate_dataset(
    t_max = t_max,
    nj = nj,
    k = k,
    mu = mu,
    lambda = lambda,
    theta = theta,
    m = m,
    beta = beta,
    A = A
  )
  
  dat
}

#dat<-Generate(condition) #block test

# =========================================================
# Analyse
# =========================================================

Analyse <- function(condition, dat, fixed_objects = NULL) {
  
  # This step is common to all scenarios
  identify_repeated_positions <- function(data, unit_id) {
    
    # Filter data by unit
    unit_data <- data %>% filter(id == unit_id)
    
    # Extract time vector
    time_vector <- unit_data$time
    
    # Compute frequency of each time value
    frequency <- table(time_vector)
    
    # Identify repeated time values
    repeated_values <- names(frequency[frequency > 1])
    
    # Identify positions of repeated values
    repeated_positions <- unlist(lapply(repeated_values, function(value) {
      which(time_vector == value)
    }))
    
    # Sort positions
    repeated_positions <- sort(repeated_positions)
    
    # Keep only odd positions
    positions <- repeated_positions[seq_along(repeated_positions) %% 2 != 0]
    
    return(positions)
  }
  
  # Compute quantities required for likelihood estimation
  # This step is also common to all units
  compute_estimation_measures <- function(data_z, positions) {
    
    # Initialize result vectors
    all_Z_increments <- c()
    all_time_increments <- c()
    all_jumps <- NULL
    all_theta <- NULL
    
    # Number of units
    m <- length(table(data_z$id))
    
    # Number of observations between maintenance actions
    nj <- positions[1] - 2
    
    # Loop through all units
    for (i in 1:m) {
      
      # Current unit subset
      subset_data <- subset(data_z, id == paste0("OBJ_", i))
      
      Z <- subset_data$Z
      t <- subset_data$time
      
      # Increments
      Z_increment <- Z[-1] - Z[-length(Z)]
      time_increment <- t[-1] - t[-length(t)]
      
      # Remove maintenance jumps
      Z_increment_clean <- Z_increment[-c(positions)]
      time_increment_clean <- time_increment[-c(positions)]
      
      # Maintenance jumps
      jump <- Z_increment[c(positions)]
      
      # Maintenance effects
      theta <- jump / (Z[c(positions)] - Z[c(positions) - (nj + 1)])
      
      # Update global vectors
      all_Z_increments <- c(all_Z_increments, Z_increment_clean)
      all_time_increments <- c(all_time_increments, time_increment_clean)
      all_jumps <- c(all_jumps, jump)
      all_theta <- c(all_theta, theta)
    }
    
    # Return results
    return(list(
      all_Z_increments = all_Z_increments,
      all_time_increments = all_time_increments,
      all_jumps = all_jumps,
      all_theta = all_theta
    ))
  }
  
  # Quantities used for parameter estimation
  
  data_z <- dat$df_z
  
  A_levels <- unique(data_z$A)
  
  all_Z_increments_global <- c()
  all_time_increments_global <- c()
  all_jumps_global <- NULL
  all_theta_global <- NULL
  A_global <- NULL
  
  for (j in 1:length(A_levels)) {
    
    A_value <- A_levels[j]
    
    # Filter data for the current acceleration level
    data_z_filtered <- data_z %>% filter(A == A_value)
    
    positions <- identify_repeated_positions(
      data_z_filtered,
      "OBJ_1"
    )
    
    estimation_measures <- compute_estimation_measures(
      data_z_filtered,
      positions
    )
    
    all_Z_increments <- estimation_measures$all_Z_increments
    all_time_increments <- estimation_measures$all_time_increments
    all_jumps <- estimation_measures$all_jumps
    all_theta <- estimation_measures$all_theta
    
    A_rep <- rep(A_value, length(all_Z_increments))
    
    # Update global vectors
    all_Z_increments_global <- c(
      all_Z_increments_global,
      all_Z_increments
    )
    
    all_time_increments_global <- c(
      all_time_increments_global,
      all_time_increments
    )
    
    all_jumps_global <- c(
      all_jumps_global,
      all_jumps
    )
    
    all_theta_global <- c(
      all_theta_global,
      all_theta
    )
    
    A_global <- c(
      A_global,
      A_rep
    )
  }
  
  
  # Data used for estimation
  estimation_data <- data.frame(
    increment = all_Z_increments_global,
    dt = all_time_increments_global,
    A = A_global
  )
  
  # Negative log-likelihood
  LogLik <- function(par, data, par.fix = kULL) {
    
    # Mean increment under accelerated degradation
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
  
  data <- estimation_data
  
  start_values <- c(1, 1, 0)
  method <- "BFGS"
  
  mod <- tryCatch(
    optim(
      fn = LogLik,
      par = start_values,
      hessian = TRUE,
      method = method,
      data = data
    ),
    error = function(e) {e}
  )
  
  e.mu <- mod$par[1]
  e.lambda <- mod$par[2]
  e.beta <- mod$par[3]
  
  mcov <- solve(mod$hessian)
  
  sd.mu <- sqrt(mcov[1,1])
  sd.lambda <- sqrt(mcov[2,2])
  sd.beta <- sqrt(mcov[3,3])
  
  IC_mu <- e.mu + c(-1,1) * qnorm(0.975) * sd.mu
  IC_lambda <- e.lambda + c(-1,1) * qnorm(0.975) * sd.lambda
  IC_beta <- e.beta + c(-1,1) * qnorm(0.975) * sd.beta
  
  results <- c(
    e.mu, sd.mu, IC_mu[1], IC_mu[2], IC_mu[2] - IC_mu[1],
    e.lambda, sd.lambda, IC_lambda[1], IC_lambda[2], IC_lambda[2] - IC_lambda[1],
    e.beta, sd.beta, IC_beta[1], IC_beta[2], IC_beta[2] - IC_beta[1]
  )
  
  results <- results %>%
    setNames(
      c(
        "mle_mu", "sd_mu", "CI_lower_mu", "CI_upper_mu", "CI_width_mu",
        "mle_lambda", "sd_lambda", "CI_lower_lambda", "CI_upper_lambda", "CI_width_lambda",
        "mle_beta", "sd_beta", "CI_lower_beta", "CI_upper_beta", "CI_width_beta"
      )
    )
  
  results
}

#results<-Analyse(condition, dat) #block test
#resultst<-rbind(results, results) #block test


# =========================================================
# Summarise (RMSE, MRE, CP)
# =========================================================

Summarise <- function(condition, results, fixed_objects = NULL) {
  
  results1 <- results
  
  mean_mu <- mean(results1[,1])
  mean_lambda <- mean(results1[,6])
  mean_beta <- mean(results1[,11])
  
  bias_mu <- mean(results1[,1] - condition$mu)
  bias_lambda <- mean(results1[,6] - condition$lambda)
  bias_beta <- mean(results1[,11] - condition$beta)
  
  RMSE_mu <- sqrt(mean((results1[,1] - condition$mu)^2))
  RMSE_lambda <- sqrt(mean((results1[,6] - condition$lambda)^2))
  RMSE_beta <- sqrt(mean((results1[,11] - condition$beta)^2))
  
  MRE_mu <- mean(results1[,1] / condition$mu)
  MRE_lambda <- mean(results1[,6] / condition$lambda)
  MRE_beta <- mean(results1[,11] / condition$beta)
  
  mean_CI_width_mu <- mean(results1[,5])
  mean_CI_width_lambda <- mean(results1[,10])
  mean_CI_width_beta <- mean(results1[,15])
  
  coverage_mu <- sum(
    ifelse(condition$mu > results1[,3] &
             condition$mu < results1[,4], 1, 0)
  )
  
  coverage_lambda <- sum(
    ifelse(condition$lambda > results1[,8] &
             condition$lambda < results1[,9], 1, 0)
  )
  
  coverage_beta <- sum(
    ifelse(condition$beta > results1[,13] &
             condition$beta < results1[,14], 1, 0)
  )
  
  n_replications <- nrow(results1)
  
  summary_results <- c(
    mean_mu = mean_mu,
    mean_lambda = mean_lambda,
    mean_beta = mean_beta,
    
    bias_mu = bias_mu,
    bias_lambda = bias_lambda,
    bias_beta = bias_beta,
    
    RMSE_mu = RMSE_mu,
    RMSE_lambda = RMSE_lambda,
    RMSE_beta = RMSE_beta,
    
    MRE_mu = MRE_mu,
    MRE_lambda = MRE_lambda,
    MRE_beta = MRE_beta,
    
    mean_CI_width_mu = mean_CI_width_mu,
    mean_CI_width_lambda = mean_CI_width_lambda,
    mean_CI_width_beta = mean_CI_width_beta,
    
    CP_mu = coverage_mu / n_replications,
    CP_lambda = coverage_lambda / n_replications,
    CP_beta = coverage_beta / n_replications
  )
  
  summary_results
}


#Summarise(condition, resultst, fixed_objects = NULL) #block test

# =========================================================
# Run Simulation
# =========================================================

final_results <- runSimulation(design=Design, replications = 5000, generate=Generate,
                           analyse=Analyse, summarise=Summarise)

final_results