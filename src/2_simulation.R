args <- commandArgs(trailingOnly = TRUE)
job_id <- if(length(args) > 0) args[1] else "0"

library(ggplot2)
library(tidyverse)
library(SuperLearner)


simulate_vaccination_study <- function(n = 500, 
                                       study_type = "RCT") {
  
  age <- rnorm(n, mean = 0, sd = 1) 
  comorbidities <- rnorm(n, mean = 0, sd = 1)
  prior_immunity <- rnorm(n, mean = 0, sd = 1)
  
  # Treatment assignment mechanism
  if (study_type == "RCT") {
    # Random assignment
    logit_ps <- rep(0, n)  # 50% probability
  } else {
    # Observational study: vaccination depends on covariates
    logit_ps <- 0.1*age + 0.01*comorbidities + 1.2*age*prior_immunity - 0.05*prior_immunity
  }
  
  ps <- plogis(logit_ps)
  A <- rbinom(n, 1, ps)
  
  # Generate potential outcomes
  # Base risk function
  
  base_risk <- plogis(
    -0.2 + 0.5*age + 0.1*comorbidities + 0.8*age*prior_immunity - 0.2*prior_immunity + 0.12 + 0.28
    
  )
  
  treatment_effect <- -0.7
  
  # Potential outcomes
  Y0 <- rbinom(n, 1, base_risk)
  Y1 <- rbinom(n, 1, plogis(qlogis(base_risk) + treatment_effect))
  
  # Observed outcome
  Y <- ifelse(A == 1, Y1, Y0)
  
  return(tibble(
    age = age, comorbidities = comorbidities, prior_immunity = prior_immunity,
    A = A, Y = Y, Y_0 = Y0, Y_1 = Y1, ps = ps
  ))
}

simulation <- simulate_vaccination_study(n = 100000, "OBS")
ATE.SOURCE.RD <- mean(simulation$Y_1) - mean(simulation$Y_0)
ATE.SOURCE.RR <- mean(simulation$Y_1) / mean(simulation$Y_0)
ATE.SOURCE.ERR <- ATE.SOURCE.RD / mean(simulation$Y_0)
ATE.SOURCE.OR <- (mean(simulation$Y_1) / (1 - mean(simulation$Y_1))) / (mean(simulation$Y_0) / (1 - mean(simulation$Y_0)))
ATE.SOURCE.NNT <- abs(1 / ATE.SOURCE.RD)
ATE.SOURCE.SR <- (1 - mean(simulation$Y_1)) / (1 - mean(simulation$Y_0))
ATE.SOURCE.VE <- 1 - ATE.SOURCE.RR
cat("Source population true RD = ", ATE.SOURCE.RD, " RR = ", ATE.SOURCE.RR, " ERR = ", ATE.SOURCE.ERR, " OR = ", ATE.SOURCE.OR, " NNT = ", ATE.SOURCE.NNT, " SR = ", ATE.SOURCE.SR, " VE = ", ATE.SOURCE.VE, "\n")

source.values <- data.frame(
  measure = c("ERR", "RD", "RR", "OR", "NNT", "SR","VE"),
  hline = c(ATE.SOURCE.ERR, ATE.SOURCE.RD, ATE.SOURCE.RR, ATE.SOURCE.OR, ATE.SOURCE.NNT, ATE.SOURCE.SR, ATE.SOURCE.VE)
)

estimate_ps <- function(train_data, test_data, correct_spec = TRUE) {
  if (correct_spec) {
    lg_fit_1 <- glm(A ~ comorbidities + age*prior_immunity, 
                    data = train_data, family = binomial())
    
    X_test <- test_data %>% dplyr::select(age, comorbidities, prior_immunity)
    ps_hat <- predict(lg_fit_1, newdata = test_data, type = "response")
    
  } else {
    ps_model_miss <- lm(A ~ age + comorbidities + prior_immunity, 
                        data = train_data)
    
    ps_hat <- predict(ps_model_miss, newdata = test_data)
  }
  
  return(ps_hat)
}

estimateSL_ps <- function(train_data, test_data, correct_spec = FALSE) {
  #SL.library <- c("SL.glm", "SL.ranger", "SL.earth", "SL.glm.interaction")
  SL.library <- c("SL.earth", "SL.glm.interaction")
  
  covars <- if (correct_spec) {
    c("age", "comorbidities", "prior_immunity")
  } else {
    c("age", "comorbidities", "prior_immunity")
  }
  
  X_train <- train_data %>% dplyr::select(all_of(covars))
  Y_train <- train_data$A
  
  sl_fit <- SuperLearner(Y = Y_train, X = X_train,
                         family = binomial(),
                         SL.library = SL.library,
                         verbose = FALSE)
  
  X_test <- test_data %>% dplyr::select(all_of(covars))
  ps_hat <- predict(sl_fit, newdata = X_test)$pred
  
  return(ps_hat)
}
estimate_outcomes <- function(train_data, test_data, correct_spec = TRUE) {
  train_data_1 <- train_data %>% filter(A == 1)
  train_data_0 <- train_data %>% filter(A == 0)
  
  if (correct_spec == TRUE) {
    # Model for A=1
    X_1 <- train_data_1 %>% dplyr::select(age, comorbidities, prior_immunity)
    Y_1 <- train_data_1$Y
    # Fit logistic regression
    lg_fit_1 <- glm(Y ~ comorbidities + age*prior_immunity, 
                    data = train_data_1, family = binomial())
    
    # Model for A=0
    X_0 <- train_data_0 %>% dplyr::select(age, comorbidities, prior_immunity) 
    Y_0 <- train_data_0$Y
    
    lg_fit_0 <- glm(Y ~ comorbidities + age*prior_immunity, 
                    data = train_data_0, family = binomial())
    # Predict on all data
    X_all <- test_data %>% dplyr::select(age, comorbidities, prior_immunity) 
    mu_1_hat <- predict(lg_fit_1, newdata = test_data, type = "response")
    mu_0_hat <- predict(lg_fit_0, newdata = test_data, type = "response")
  } else {
    # Misspecified models - linear only, missing interactions
    mu_1_model_miss <- lm(Y ~ age + comorbidities + prior_immunity, data = train_data_1)
    mu_0_model_miss <- lm(Y ~ age + comorbidities + prior_immunity, data = train_data_0)
    
    mu_1_hat <- predict(mu_1_model_miss, newdata = test_data)
    mu_0_hat <- predict(mu_0_model_miss, newdata = test_data)
  }
  
  return(list(mu_1 = mu_1_hat, mu_0 = mu_0_hat))
}

estimateSL_outcomes <- function(train_data, test_data, correct_spec = TRUE) {
  train_data_1 <- train_data %>% filter(A == 1)
  train_data_0 <- train_data %>% filter(A == 0)
  
  SL.library <- c("SL.earth", "SL.glm.interaction")
  #SL.library <- c("SL.glm", "SL.ranger", "SL.earth", "SL.glm.interaction")
  
  if (correct_spec == TRUE) {
    covars <- c("age", "comorbidities", "prior_immunity")
  } else {
    covars <- c("age", "comorbidities", "prior_immunity")
  }
  
  # Model for A=1
  X_1 <- train_data_1 %>% dplyr::select(all_of(covars))
  Y_1 <- train_data_1$Y
  sl_fit_1 <- SuperLearner(Y = Y_1, X = X_1,
                           family = binomial(),
                           SL.library = SL.library,
                           verbose = FALSE)
  
  # Model for A=0
  X_0 <- train_data_0 %>% dplyr::select(all_of(covars))
  Y_0 <- train_data_0$Y
  sl_fit_0 <- SuperLearner(Y = Y_0, X = X_0,
                           family = binomial(),
                           SL.library = SL.library,
                           verbose = FALSE)
  
  # Predict on test data
  X_test <- test_data %>% dplyr::select(all_of(covars))
  mu_1_hat <- predict(sl_fit_1, newdata = X_test)$pred
  mu_0_hat <- predict(sl_fit_0, newdata = X_test)$pred
  
  return(list(mu_1 = as.vector(mu_1_hat), mu_0 = as.vector(mu_0_hat)))
}

compute_aipw_cf <- function(data, ps_correct = TRUE, outcome_correct = TRUE, n_folds = 5) {
  n <- nrow(data)
  fold_ids <- sample(rep(1:n_folds, length.out = n))
  aipw_1_n <- numeric(n)
  aipw_0_n <- numeric(n)
  
  for (k in 1:n_folds) {
    # Split data
    train_idx <- fold_ids != k
    test_idx <- fold_ids == k
    
    train_data <- data[train_idx, ]
    test_data <- data[test_idx, ]
    
    # Estimate propensity scores
    ps_hat <- estimate_ps(train_data, test_data, correct_spec = ps_correct)
    
    # Estimate outcome models
    mu_hats <- estimate_outcomes(train_data, test_data, correct_spec = outcome_correct)
    mu_1_hat <- mu_hats$mu_1
    mu_0_hat <- mu_hats$mu_0
    
    ps_hat <- pmax(ps_hat, 0.01)
    ps_hat <- pmin(ps_hat, 0.99)
    
    aipw_1 <- mu_1_hat + (test_data$A / ps_hat) * (test_data$Y - mu_1_hat)
    aipw_0 <- mu_0_hat + ((1 - test_data$A) / (1 - ps_hat)) * (test_data$Y - mu_0_hat)
    
    aipw_1_n[test_idx] <- aipw_1
    aipw_0_n[test_idx] <- aipw_0
    
  }
  
  psi_1 <- mean(aipw_1_n)  # E[Y(1)]
  psi_0 <- mean(aipw_0_n)
  
  RD <- psi_1 - psi_0
  RR <- psi_1 / psi_0
  ERR <- (psi_1 - psi_0) / psi_0
  OR <- (psi_1 / (1 - psi_1)) /(psi_0 / (1 - psi_0))
  NNT <- abs(1 / (psi_1 - psi_0))
  SR <- (1 - psi_1) / (1 - psi_0)
  VE <- 1 - psi_1 / psi_0
  
  results_df <- data.frame(
    measure = c("RD", "ERR","RR",  "OR", "NNT", "SR", "VE"),
    estimate = c(RD, RR, ERR, OR, NNT, SR, VE)
  )
  
  IF_psi1 <- aipw_1_n - psi_1
  
  IF_psi0 <- aipw_0_n - psi_0
  
  # Influence functions for RD
  IF_RD <- IF_psi1 - IF_psi0
  
  # Influence functions for log RR
  IF_logRR <- (1/psi_1) * IF_psi1 - 
    (1/psi_0) * IF_psi0
  
  # Influence functions for log ERR
  IF_logERR <- (1/(psi_1 - psi_0)) * IF_psi1 + 
    (-1/(psi_1 - psi_0) - 1/psi_0) * IF_psi0
  
  # Influence functions for log OR
  IF_logOR <- (1/psi_1 + 1/(1 - psi_1)) * IF_psi1 + 
    (-1/(1 - psi_0) - 1/psi_0) * IF_psi0
  
  # Influence functions for log NNT
  IF_logNNT <- (1/(psi_0 - psi_1)) * IF_psi1 + 
    (-1/(psi_0 - psi_1)) * IF_psi0
  
  # Influence functions for log SR
  IF_logSR <- (-1/(1 - psi_1)) * IF_psi1 + 
    (1/(1 - psi_0)) * IF_psi0
  
  # Influence functions for log VE
  IF_logVE <- (-1/(psi_0 - psi_1)) * IF_psi1 + 
    (1/(psi_0 - psi_1) - 1/psi_0) * IF_psi0
  
  IF_log_list <- list(IF_logRR, IF_logOR, IF_logNNT, IF_logSR)
  se_log <- sapply(IF_log_list, function(x) sqrt(var(x) / n))
  se_RD <- sqrt(var(IF_RD) / n)
  z <- qnorm(0.975)
  RD_lower <- RD - z * se_RD
  RD_upper <- RD + z * se_RD
  
  
  log_estimates <- log(c(RR, OR, NNT, SR))
  log_lower <- log_estimates - z * se_log
  log_upper <- log_estimates + z * se_log
  
  # ERR = RR - 1: derive from RR
  ERR_lower <- exp(log_lower[1]) - 1
  ERR_upper <- exp(log_upper[1]) - 1
  
  # VE = 1 - RR: derive from RR (bounds flip)
  VE_lower <- 1 - exp(log_upper[1])
  VE_upper <- 1 - exp(log_lower[1])
  
  lower <- c(RD_lower, ERR_lower, exp(log_lower[1]), exp(log_lower[2]), 
             exp(log_lower[3]), exp(log_lower[4]), VE_lower)
  upper <- c(RD_upper, ERR_upper, exp(log_upper[1]), exp(log_upper[2]), 
             exp(log_upper[3]), exp(log_upper[4]), VE_upper)
  
  results_df <- data.frame(
    measure  = c("RD", "ERR", "RR", "OR", "NNT", "SR", "VE"),
    estimates = c(RD, ERR, RR, OR, NNT, SR, VE),
    lower,
    upper
  )
  
  
  return(results_df)
}

computeSL_aipw_cf <- function(data, ps_correct = TRUE, outcome_correct = TRUE, n_folds = 5) {
  n <- nrow(data)
  fold_ids <- sample(rep(1:n_folds, length.out = n))
  aipw_1_n <- numeric(n)
  aipw_0_n <- numeric(n)
  
  for (k in 1:n_folds) {
    # Split data
    train_idx <- fold_ids != k
    test_idx <- fold_ids == k
    
    train_data <- data[train_idx, ]
    test_data <- data[test_idx, ]
    
    # Estimate propensity scores
    ps_hat <- estimateSL_ps(train_data, test_data, correct_spec = ps_correct)
    
    # Estimate outcome models
    mu_hats <- estimateSL_outcomes(train_data, test_data, correct_spec = outcome_correct)
    mu_1_hat <- mu_hats$mu_1
    mu_0_hat <- mu_hats$mu_0
    
    ps_hat <- pmax(ps_hat, 0.01)
    ps_hat <- pmin(ps_hat, 0.99)
    
    aipw_1 <- mu_1_hat + (test_data$A / ps_hat) * (test_data$Y - mu_1_hat)
    aipw_0 <- mu_0_hat + ((1 - test_data$A) / (1 - ps_hat)) * (test_data$Y - mu_0_hat)
    
    aipw_1_n[test_idx] <- aipw_1
    aipw_0_n[test_idx] <- aipw_0
    
  }
  
  psi_1 <- mean(aipw_1_n)  # E[Y(1)]
  psi_0 <- mean(aipw_0_n)
  
  RD <- psi_1 - psi_0
  RR <- psi_1 / psi_0
  ERR <- (psi_1 - psi_0) / psi_0
  OR <- (psi_1 / (1 - psi_1)) /(psi_0 / (1 - psi_0))
  NNT <- abs(1 / (psi_1 - psi_0))
  SR <- (1 - psi_1) / (1 - psi_0)
  VE <- 1 - psi_1 / psi_0
  
  results_df <- data.frame(
    measure = c("RD", "RR", "ERR", "OR", "NNT", "SR", "VE"),
    estimate = c(RD, RR, ERR, OR, NNT, SR, VE)
  )
  
  IF_psi1 <- aipw_1_n - psi_1
  
  IF_psi0 <- aipw_0_n - psi_0
  
  # Influence functions for RD
  IF_RD <- IF_psi1 - IF_psi0
  
  
  # Influence functions for log RR
  IF_logRR <- (1/psi_1) * IF_psi1 - 
    (1/psi_0) * IF_psi0
  
  # Influence functions for log ERR
  IF_logERR <- (1/(psi_1 - psi_0)) * IF_psi1 + 
    (-1/(psi_1 - psi_0) - 1/psi_0) * IF_psi0
  
  # Influence functions for log OR
  IF_logOR <- (1/psi_1 + 1/(1 - psi_1)) * IF_psi1 + 
    (-1/(1 - psi_0) - 1/psi_0) * IF_psi0
  
  # Influence functions for log NNT
  IF_logNNT <- (1/(psi_0 - psi_1)) * IF_psi1 + 
    (-1/(psi_0 - psi_1)) * IF_psi0
  
  # Influence functions for log SR
  IF_logSR <- (-1/(1 - psi_1)) * IF_psi1 + 
    (1/(1 - psi_0)) * IF_psi0
  
  # Influence functions for log VE
  IF_logVE <- (-1/(psi_0 - psi_1)) * IF_psi1 + 
    (1/(psi_0 - psi_1) - 1/psi_0) * IF_psi0
  
  IF_log_list <- list(IF_logRR, IF_logOR, IF_logNNT, IF_logSR)
  se_log <- sapply(IF_log_list, function(x) sqrt(var(x) / n))
  se_RD <- sqrt(var(IF_RD) / n)
  z <- qnorm(0.975)
  RD_lower <- RD - z * se_RD
  RD_upper <- RD + z * se_RD
  
  log_estimates <- log(c(RR, OR, NNT, SR))
  log_lower <- log_estimates - z * se_log
  log_upper <- log_estimates + z * se_log
  
  # ERR = RR - 1: derive from RR
  ERR_lower <- exp(log_lower[1]) - 1
  ERR_upper <- exp(log_upper[1]) - 1
  
  # VE = 1 - RR: derive from RR (bounds flip)
  VE_lower <- 1 - exp(log_upper[1])
  VE_upper <- 1 - exp(log_lower[1])
  
  lower <- c(RD_lower, ERR_lower, exp(log_lower[1]), exp(log_lower[2]), 
             exp(log_lower[3]), exp(log_lower[4]), VE_lower)
  upper <- c(RD_upper, ERR_upper, exp(log_upper[1]), exp(log_upper[2]), 
             exp(log_upper[3]), exp(log_upper[4]), VE_upper)
  
  results_df <- data.frame(
    measure  = c("RD", "ERR", "RR", "OR", "NNT", "SR", "VE"),
    estimates = c(RD, ERR, RR, OR, NNT, SR, VE),
    lower,
    upper
  )
  
  
  return(results_df)
}
#b=200

compute_aipw_bootstrap <- function(data, ps_correct = TRUE, outcome_correct = TRUE, n_boot = 200) {
  n <- nrow(data)
  
  # Original estimate
  orig_est <- compute_aipw_cf(data, ps_correct, outcome_correct)
  
  # Bootstrap estimates
  boot_estimates <- matrix(NA, nrow = n_boot, ncol = 7)
  colnames(boot_estimates) <- c("RD", "ERR", "RR", "OR", "NNT", "SR", "VE")
  
  for (b in 1:n_boot) {
    boot_indices <- sample(1:n, size = n, replace = TRUE)
    boot_data <- data[boot_indices, ]
    boot_result <- compute_aipw_cf(boot_data, ps_correct, outcome_correct)
    boot_estimates[b, ] <- boot_result$estimates
  }
  
  log_measures <- c("RR", "OR", "NNT", "SR")
  log_boot <- log(boot_estimates[, log_measures])
  log_lower <- apply(log_boot, 2, quantile, probs = 0.025, na.rm = TRUE)
  log_upper <- apply(log_boot, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  # RD: directly on original scale
  RD_lower <- quantile(boot_estimates[, "RD"], 0.025, na.rm = TRUE)
  RD_upper <- quantile(boot_estimates[, "RD"], 0.975, na.rm = TRUE)
  
  # ERR = RR - 1: derive from RR
  ERR_lower <- exp(log_lower["RR"]) - 1
  ERR_upper <- exp(log_upper["RR"]) - 1
  
  # VE = 1 - RR: derive from RR (bounds flip)
  VE_lower <- 1 - exp(log_upper["RR"])
  VE_upper <- 1 - exp(log_lower["RR"])
  
  boot_lower <- c(RD_lower, ERR_lower, exp(log_lower["RR"]), exp(log_lower["OR"]), 
                  exp(log_lower["NNT"]), exp(log_lower["SR"]), VE_lower)
  boot_upper <- c(RD_upper, ERR_upper, exp(log_upper["RR"]), exp(log_upper["OR"]), 
                  exp(log_upper["NNT"]), exp(log_upper["SR"]), VE_upper)
  
  results_df <- data.frame(
    measure  = c("RD", "ERR", "RR", "OR", "NNT", "SR", "VE"),
    estimates = orig_est$estimates,
    lower    = boot_lower,
    upper    = boot_upper
  )
  rownames(results_df) <- NULL
  
  return(results_df)
}

n_sample <- 5000
n_sim <- 1

# Initialize results storage
results <- data.frame(
  measure = character(),
  setting = character(),
  estimates = numeric(),
  lower = numeric(),
  upper = numeric()
)

# Generate data
sim_data <- simulate_vaccination_study(n = n_sample, "OBS")

# Setting 1: Both correctly specified
aipw_1 <- compute_aipw_cf(sim_data, ps_correct = TRUE, outcome_correct = TRUE)

# Setting 2: PS correct, outcome misspecified
aipw_2 <- compute_aipw_bootstrap(sim_data, ps_correct = TRUE, outcome_correct = FALSE, n_boot = 500)

# Setting 3: PS misspecified, outcome correct
aipw_3 <- compute_aipw_bootstrap(sim_data, ps_correct = FALSE, outcome_correct = TRUE, n_boot = 500)

aipw_4 <- compute_aipw_bootstrap(sim_data, ps_correct = FALSE, outcome_correct = FALSE, n_boot = 500)

aipw_5 <- computeSL_aipw_cf(sim_data, ps_correct = FALSE, outcome_correct = FALSE)

# Store results
settings <- c("Both Correct", "PS Correct, Outcome Wrong", 
              "PS Wrong, Outcome Correct", "PS Wrong, Outcome Wrong",
              "SuperLearner with all variables")
aipw_list <- list(aipw_1, aipw_2, aipw_3, aipw_4, aipw_5)

for (j in seq_along(settings)) {
  df <- aipw_list[[j]]
  df$setting <- settings[j]
  results <- rbind(results, df)
}

output_filename <- paste0("simulation_results_", job_id, ".csv")
write.csv(results, output_filename, row.names = FALSE)
print(paste("Results saved to:", output_filename))
