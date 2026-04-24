# =========================================================
# Inverse Gaussian Degradation with Exponential Acceleration
# Toy Example
# =========================================================

# Packages
library(statmod)
library(ggplot2)

# =========================================================
# Function: Generate one degradation trajectory
# =========================================================
generate_trajectory <- function(t_max, nj, k, mu, lambda, beta, s_value, s_label, theta, id) {
  
  # Time grid
  N <- nj * (k + 1)
  n_points <- N + k + 2
  time_x <- seq(0, t_max, length.out = n_points)
  dt <- t_max / (n_points - 1)
  
  # Inverse Gaussian increments with exponential acceleration
  increments <- rinvgauss(
    n_points - 1,
    mean = mu * dt * exp(beta * s_value),
    shape = (mu * dt * exp(beta * s_value))^2 * lambda
  )
  
  # Latent degradation process X(t)
  X <- c(0, cumsum(increments))
  
  # Maintenance times
  t_change <- seq(0, t_max, length.out = k + 2)
  t_MP <- t_change[-c(1, length(t_change))]
  
  # Observed process Z(t) includes maintenance effects
  time_z <- sort(c(time_x, t_MP))
  Z <- numeric(length(time_z))
  
  j <- length(Z) / (k + 1)
  l <- (length(X) - 1) / (k + 1)
  
  # Apply imperfect maintenance (theta)
  for (i in 1:k) {
    Z[1:j] <- X[1:(l + 1)]
    
    reduction <- sum(sapply(1:i, function(p) {
      theta[p] * (X[p * l + 1] - X[(p - 1) * l + 1])
    }))
    
    Z[(i * j + 1):((i + 1) * j)] <- 
      X[(i * l + 1):((i + 1) * l + 1)] - reduction
  }
  
  # Output data frames
  df_x <- data.frame(
    id = id,
    time = time_x,
    X = X,
    A = s_label
  )
  
  df_z <- data.frame(
    id = id,
    time = time_z,
    Z = Z,
    A = s_label
  )
  
  return(list(df_x = df_x, df_z = df_z, t_MP = t_MP))
}

# =========================================================
# Function: Generate full dataset across stress levels
# =========================================================
generate_dataset <- function(n_units = 5,
                             A_values = c(0, 0.15, 0.3),
                             A_labels = c("A1", "A2", "A3"),
                             t_max = 20,
                             nj = 2,
                             k = 4,
                             mu = 3,
                             lambda = 2,
                             beta = 2,
                             theta = c(0.9, 0.4, 0.6, 0.9)) {
  
  all_x <- list()
  all_z <- list()
  
  for (j in seq_along(A_values)) {
    for (i in 1:n_units) {
      
      sim <- generate_trajectory(
        t_max, nj, k, mu, lambda, beta,
        s_value = A_values[j],
        s_label = A_labels[j],
        theta = theta,
        id = i
      )
      
      all_x[[length(all_x) + 1]] <- sim$df_x
      all_z[[length(all_z) + 1]] <- sim$df_z
    }
  }
  
  df_x <- do.call(rbind, all_x)
  df_z <- do.call(rbind, all_z)
  
  return(list(df_x = df_x, df_z = df_z, t_MP = sim$t_MP))
}

# =========================================================
# Function: Plot Z(t) (with maintenance)
# =========================================================
plot_degradation <- function(data) {
  ggplot(data$df_z, aes(x = time, y = Z, 
                        group = interaction(id, A),
                        color = A)) +
    geom_line(alpha = 0.7) +
    geom_point(size = 1) +
    geom_vline(xintercept = data$t_MP,
               linetype = "dashed", color = "gray") +
    labs(
      x = "Time",
      y = "Z(t)",
      color = "Acceleration level"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
}

# =========================================================
# Example
# =========================================================
set.seed(123)

data_sim <- generate_dataset()

# Visualize Z(t)
plot_degradation(data_sim)

# =========================================================
# Combined plots: X(t), Z(t), and both
# =========================================================
library(patchwork)

# Z(t): degradation with maintenance
p_Z <- ggplot(data_sim$df_z,
              aes(time, Z, group = interaction(id, A), color = A)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1) +
  geom_vline(xintercept = data_sim$t_MP,
             linetype = "dashed", color = "gray") +
  labs(x = "Time", y = "Z(t)", color = "Acceleration level") +
  theme_bw() +
  theme(legend.position = "none")

# X(t): degradation without maintenance
p_X <- ggplot(data_sim$df_x,
              aes(time, X, group = interaction(id, A), color = A)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1) +
  geom_vline(xintercept = data_sim$t_MP,
             linetype = "dashed", color = "gray") +
  labs(x = "Time", y = "X(t)", color = "Acceleration level") +
  theme_bw() +
  theme(legend.position = "none")

# Combined plot: X(t) and Z(t)
p_both <- ggplot() +
  geom_line(data = data_sim$df_z,
            aes(time, Z, group = interaction(id, A), color = A),
            alpha = 0.8) +
  geom_line(data = data_sim$df_x,
            aes(time, X, group = interaction(id, A), color = A),
            alpha = 0.25) +
  geom_vline(xintercept = data_sim$t_MP,
             linetype = "dashed", color = "gray") +
  labs(x = "Time", y = "Degradation", color = "Acceleration level") +
  theme_bw() +
  theme(legend.position = "bottom")

# Combine panels vertically
final_plot <- p_Z / p_X / p_both

# Display final figure
final_plot