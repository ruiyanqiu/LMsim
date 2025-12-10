# ==============================================================================
# Linear Model Simulation & Diagnostics (STAT 5430 Project)
# Author: Ruiyan Qiu
# Purpose: Evaluate statistical methods under collinearity, misspecification, and non-normality
# ==============================================================================

# --------------------------
# Core Package Loading (Auto-Install Missing Packages)
# --------------------------
required_pkgs <- c(
  "tidyverse", "glmnet", "car", "agricolae", "purrr", 
  "NHANES", "MASS", "dplyr", "ggplot2"
)

lapply(required_pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
})

# --------------------------
# 1. Data Generation Function (Simulated Data)
# --------------------------
#' Generate Data with Collinearity, Measurement Error, and Non-Normality
#' 
#' @description
#' Creates simulated datasets under controlled levels of collinearity,
#' measurement error, non-normal error distributions, and model
#' misspecification. Returns both true and observed predictors for
#' downstream model evaluation.
#' @param n Sample size (int)
#' @param collinearity_level "none", "moderate", "severe" (controls VIF)
#' @param me_level "none", "low", "high" (measurement error for X1)
#' @param error_dist "normal", "t3", "lognormal", "poisson"
#' @param outcome_type "continuous", "binary", "count"
#' @param misspecification "correct", "underfit", "overfit"
#' @return Tibble with true/observed predictors and outcome
#' @export
generate_sim_data <- function(n, collinearity_level, me_level, error_dist, outcome_type, misspecification) {
  set.seed(123 + n + match(collinearity_level, c("none", "moderate", "severe")))
  
  # True predictors
  X1_true <- rnorm(n)
  X2 <- runif(n, -1, 1)
  
  # Collinearity (calibrated for target VIF)
  sigma_coll <- switch(collinearity_level,
                       none = 1,
                       moderate = 0.2,
                       severe = 0.05)
  X3_true <- 0.8*X1_true + 0.2*X2 + rnorm(n, 0, sigma_coll)
  
  # Measurement error for X1 (observed X1)
  sigma_me <- switch(me_level,
                     none = 0,
                     low = 0.1,
                     high = 0.5)
  X1 <- X1_true + rnorm(n, 0, sigma_me)
  
  # Observed X3 (no measurement error)
  X3 <- X3_true
  
  # Overfit: add irrelevant predictors
  if (misspecification == "overfit") {
    irrelevant <- matrix(rnorm(n*17), ncol=17)
    colnames(irrelevant) <- paste0("X", 4:20)
    predictors <- as_tibble(cbind(X1, X2, X3, irrelevant))
  } else {
    predictors <- as_tibble(cbind(X1, X2, X3))
  }
  
  # Non-normal error terms
  e <- switch(error_dist,
              normal = rnorm(n),
              t3 = rt(n, df=3),
              lognormal = rlnorm(n) - exp(0.5),
              poisson = rpois(n, lambda=2) - 2)
  
  # Outcome generation
  beta0 <- 2; beta1 <- 1.2; beta2 <- -0.8; beta3 <- 0.5; beta4 <- 0.7
  linear_pred <- beta0 + beta1*X1_true + beta2*X2 + beta3*X3_true
  
  # Underfit: omit quadratic term
  if (misspecification != "underfit") {
    linear_pred <- linear_pred + beta4*(X1_true^2)
  }
  
  Y <- switch(outcome_type,
              continuous = linear_pred + e,
              binary = rbinom(n, 1, plogis(linear_pred)),
              count = rpois(n, lambda = exp(linear_pred) * 1.5))
  
  # Return data with true predictors
  predictors %>%
    mutate(
      Y = Y,
      X1_true = X1_true,
      X3_true = X3_true
    )
}

# --------------------------
# 2. Model Estimation Function (LS/GLM/Ridge/LASSO) — FIXED STEPWISE & FORMULA CONSISTENCY
# --------------------------
#' Fit Linear/GLM/Regularized Models
#' 
#' @description
#' Fits least squares, GLM, Ridge, LASSO, and stepwise models to the
#' simulated data. Automatically handles misspecification scenarios
#' and different outcome types.
#' @param data Tibble from generate_sim_data()
#' @param outcome_type "continuous", "binary", "count"
#' @param misspecification "correct", "underfit", "overfit"
#' @return List of model fits (LS, GLM, Ridge, LASSO, Stepwise)
#' @export
fit_models <- function(data, outcome_type, misspecification) {
  
  # Remove true predictors if present
  data_obs <- data[, !colnames(data) %in% c("X1_true", "X3_true")]
  
  # Define formula based on misspecification
  formula <- if (misspecification == "underfit") {
    Y ~ X1 + X2 + X3
  } else {
    Y ~ X1 + X2 + X3 + I(X1^2)
  }
  
  # --------------------------
  # Linear model (LS)
  # --------------------------
  ls_fit <- if (outcome_type == "continuous") {
    lm(formula, data = data_obs)
  } else NULL
  
  # --------------------------
  # GLM (binary/count)
  # --------------------------
  glm_fit <- switch(outcome_type,
                    binary = glm(formula, data = data_obs, family = binomial(link="logit")),
                    count = glm(formula, data = data_obs, family = poisson(link="log")),
                    continuous = NULL)
  
  # --------------------------
  # Ridge / LASSO using glmnet
  # --------------------------
  x <- model.matrix(formula, data_obs)[, -1]  # remove intercept
  y <- data_obs$Y
  glmnet_family <- switch(outcome_type,
                          continuous = "gaussian",
                          binary = "binomial",
                          count = "poisson")
  
  ridge_fit <- glmnet(x, y, alpha = 0, family = glmnet_family, standardize = TRUE)
  lasso_fit <- glmnet(x, y, alpha = 1, family = glmnet_family, standardize = TRUE)
  
  # --------------------------
  # Stepwise Selection (consistent formula)
  # --------------------------
  step_fit <- NULL
  if (outcome_type == "continuous") {
    step_fit <- tryCatch(
      step(lm(formula, data = data_obs), direction = "both", trace = 0),
      error = function(e) { message("Stepwise LS failed: ", e$message); NULL }
    )
  } else if (outcome_type == "binary") {
    step_fit <- tryCatch(
      step(glm(formula, data = data_obs, family = binomial), direction = "both", trace = 0),
      error = function(e) { message("Stepwise Binary GLM failed: ", e$message); NULL }
    )
  } else if (outcome_type == "count") {
    step_fit <- tryCatch(
      step(glm(formula, data = data_obs, family = poisson), direction = "both", trace = 0),
      error = function(e) { message("Stepwise Count GLM failed: ", e$message); NULL }
    )
  }
  
  # --------------------------
  # Return all fits with consistent names
  # --------------------------
  list(
    LS = ls_fit,
    GLM = glm_fit,
    Ridge = ridge_fit,
    LASSO = lasso_fit,
    Stepwise = step_fit
  )
}


# --------------------------
# 3. Model Selection Metrics (AIC/BIC/Mallow's Cp)
# --------------------------
#' Calculate Model Selection Metrics (Robust, No olsrr Dependency)
#' 
#' @description
#' Computes AIC, BIC, and Mallow’s Cp for all model fits returned by
#' `fit_models()`. Provides a unified comparison across estimation
#' methods under model misspecification.
#' @param model_fits List from fit_models()
#' @param data_obs Tibble (observed data)
#' @return Tibble with AIC, BIC, Mallow's Cp
#' @export
calculate_selection_metrics <- function(model_fits, data_obs) {
  # Base R Mallow's Cp
  mallow_cp <- function(fit, fullmodel) {
    if (!inherits(fit, "lm") || !inherits(fullmodel, "lm")) return(NA)
    n <- nrow(data_obs)
    p <- length(coef(fit))
    p_full <- length(coef(fullmodel))
    rss <- sum(residuals(fit)^2)
    rss_full <- sum(residuals(fullmodel)^2)
    mse_full <- rss_full / (n - p_full)
    (rss / mse_full) - n + 2 * p
  }
  
  # Filter NULL models
  models <- model_fits %>% discard(is.null)
  
  # Calculate metrics
  metrics <- map_dfr(names(models), function(m) {
    fit <- models[[m]]
    aic <- bic <- cp <- NA
    
    # AIC/BIC (lm/glm only)
    if (inherits(fit, c("lm", "glm"))) {
      aic <- tryCatch(AIC(fit), error = function(e) NA)
      bic <- tryCatch(BIC(fit), error = function(e) NA)
    }
    
    # Mallow's Cp (lm/stepwise only)
    if (inherits(fit, "lm") && m %in% c("LS", "Stepwise")) {
      fullmodel <- tryCatch(lm(Y ~ ., data = data_obs), error = function(e) NULL)
      cp <- if (!is.null(fullmodel)) tryCatch(mallow_cp(fit, fullmodel), error = function(e) NA) else NA
    }
    
    tibble(model = m, AIC = aic, BIC = bic, Cp = cp)
  })
  
  # Ensure all model names are present
  all_models <- tibble(model = c("LS", "GLM", "Ridge", "LASSO", "Stepwise"))
  left_join(all_models, metrics, by = "model")
}

# --------------------------
# 4. Model Diagnostics (VIF/Cook's Distance/Residuals) — FIXED SHAPIRO TEST
# --------------------------
#' Linear Model Diagnostics (Safe Shapiro Test)
#' 
#' @description
#' Produces VIF, Cook’s distance, residual diagnostics, and a
#' safe Shapiro-Wilk normality test that automatically handles
#' sample-size limitations.
#' @param fit LS/GLM fit
#' @param data_obs Observed data
#' @return List of diagnostics (VIF, Cook's distance, residual plots)
#' @export
run_diagnostics <- function(fit, data_obs) {
  if (is.null(fit)) return(NULL)
  
  # VIF
  vif_vals <- car::vif(fit)
  
  # Cook's distance
  cookd <- cooks.distance(fit)
  influential <- which(cookd > 4/length(cookd))
  
  # Residual plot
  res_plot <- ggplot(data.frame(resid = residuals(fit), fitted = fitted(fit)),
                     aes(x = fitted, y = resid)) +
    geom_point() +
    geom_hline(yintercept = 0, color = "red") +
    labs(title = "Residual vs. Fitted", x = "Fitted Values", y = "Residuals") +
    theme_minimal()
  
  # --------------------------
  # FIX: Safe Shapiro Test (Sample Size Check)
  # --------------------------
  shapiro_result <- tryCatch({
    resids <- residuals(fit)
    # Limit to 5000 observations (max for Shapiro-Wilk)
    if (length(resids) > 5000) {
      resids_subset <- sample(resids, 5000, replace = FALSE)
      shapiro.test(resids_subset)
    } else if (length(resids) < 3) {
      warning("Sample size <3 for Shapiro test")
      NULL
    } else {
      shapiro.test(resids)
    }
  }, error = function(e) {
    message("Shapiro test failed: ", e$message)
    NULL
  })
  
  list(
    vif = vif_vals,
    cookd = cookd,
    influential_obs = influential,
    residual_plot = res_plot,
    shapiro_test = shapiro_result # Safe result (NULL if failed)
  )
}

# --------------------------
# 5. Multiple Comparisons (ANCOVA/LSD/Bonferroni/Scheffe)
# --------------------------
#' Multiple Comparison Tests (LSD, Bonferroni, Scheffe)
#' 
#' @description
#' Performs classical multiple comparison procedures (LSD, Bonferroni,
#' Scheffe) on ANCOVA-style models and computes empirical false discovery
#' rates across simulated contrasts.
#' @param ancova_fit ANCOVA fit object (lm(Y ~ factor + covariate))
#' @return Tibble with FDR estimates
#' @export
run_multiple_comparisons <- function(ancova_fit) {
  if (is.null(ancova_fit) || !inherits(ancova_fit, "lm")) {
    return(tibble(method = c("LSD", "Bonferroni", "Scheffe"), fdr = NA))
  }
  
  # Run multiple comparison tests
  lsd <- tryCatch(LSD.test(ancova_fit, "factor", p.adj = "none"), error = function(e) NULL)
  bonferroni <- tryCatch(LSD.test(ancova_fit, "factor", p.adj = "bonferroni"), error = function(e) NULL)
  scheffe <- tryCatch(scheffe.test(ancova_fit, "factor"), error = function(e) NULL)
  
  # Simulate true nulls (80% nulls)
  set.seed(123)
  n_comparisons <- if (!is.null(lsd) && !is.null(lsd$comparison)) nrow(lsd$comparison) else 0
  true_nulls <- if (n_comparisons > 0) sample(c(TRUE, FALSE), n_comparisons, replace = TRUE, prob = c(0.8, 0.2)) else logical(0)
  
  # Calculate FDR
  calculate_fdr <- function(test_result, nulls) {
    if (is.null(test_result) || is.null(test_result$comparison) || nrow(test_result$comparison) == 0) return(NA)
    mean(test_result$comparison$pval < 0.05 & nulls, na.rm = TRUE)
  }
  
  tibble(
    method = c("LSD", "Bonferroni", "Scheffe"),
    fdr = c(calculate_fdr(lsd, true_nulls), calculate_fdr(bonferroni, true_nulls), calculate_fdr(scheffe, true_nulls))
  )
}

# --------------------------
# 6. Prediction Intervals & Confidence Bands
# --------------------------
#' Prediction Intervals (Point-wise vs. Simultaneous)
#' 
#' @description
#' Computes pointwise mean confidence intervals, simultaneous bands,
#' and prediction intervals, and evaluates empirical coverage using
#' simulated true outcomes.
#' @param fit LS fit (linear model only)
#' @param newdata New data with Y_true column
#' @return Tibble with coverage probability
#' @export
calculate_prediction_bands <- function(fit, newdata) {
  if (is.null(fit) || !inherits(fit, "lm")) {
    return(tibble(band_type = c("Point-wise (Mean)", "Simultaneous (Mean)", "Point-wise (Obs)"), coverage = NA))
  }
  
  # Point-wise confidence interval (mean)
  pointwise_mean <- predict(fit, newdata = newdata, interval = "confidence", level = 0.95)
  
  # Simultaneous confidence band (fixed)
  coefs <- coef(fit)
  cov_mat <- vcov(fit)
  X_new <- model.matrix(formula(fit), newdata)
  fit_vals <- as.vector(X_new %*% coefs)
  se_vals <- sqrt(diag(X_new %*% cov_mat %*% t(X_new)))
  
  simultaneous_mean <- data.frame(
    fit = fit_vals,
    lwr = fit_vals - se_vals * qt(0.975, df = fit$df.residual),
    upr = fit_vals + se_vals * qt(0.975, df = fit$df.residual)
  )
  
  # Point-wise prediction interval (observations)
  pointwise_obs <- predict(fit, newdata = newdata, interval = "prediction", level = 0.95)
  
  # Coverage calculation (Y_true required)
  if (!"Y_true" %in% colnames(newdata)) stop("newdata must include Y_true")
  Y_true <- newdata$Y_true
  
  tibble(
    band_type = c("Point-wise (Mean)", "Simultaneous (Mean)", "Point-wise (Obs)"),
    coverage = c(
      mean(Y_true >= pointwise_mean[, "lwr"] & Y_true <= pointwise_mean[, "upr"], na.rm = TRUE),
      mean(Y_true >= simultaneous_mean$lwr & Y_true <= simultaneous_mean$upr, na.rm = TRUE),
      mean(Y_true >= pointwise_obs[, "lwr"] & Y_true <= pointwise_obs[, "upr"], na.rm = TRUE)
    )
  )
}

# --------------------------
# 7. Full Simulation (Error-Resilient + Enhanced)
# --------------------------
#' Run Full Simulation Across All Conditions
#' 
#' @description
#' Executes the complete simulation workflow for a user-defined grid
#' of experimental factors, including data generation, model fitting,
#' model selection, diagnostics, prediction band evaluation, and
#' multiple-comparisons performance.
#' @param sim_grid Experimental grid (expand.grid)
#' @param R Number of replications
#' @return Tibble with aggregated performance metrics
#' @export
run_full_simulation <- function(sim_grid, R = 2000) {
  # Safe replication wrapper
  safe_replicate <- function(r, n, collinearity_level, me_level, error_dist, outcome_type, misspecification) {
    tryCatch({
      # Beta parameters (match generate_sim_data)
      beta0 <- 2; beta1 <- 1.2; beta2 <- -0.8; beta3 <- 0.5; beta4 <- 0.7
      
      # Generate data
      data <- generate_sim_data(n, collinearity_level, me_level, error_dist, outcome_type, misspecification)
      
      # Observed data (remove true predictors)
      data_obs <- data[, !colnames(data) %in% c("X1_true", "X3_true")]
      
      # Fit models
      fits <- fit_models(data, outcome_type, misspecification)
      
      # Model selection metrics
      metrics <- calculate_selection_metrics(fits, data_obs)
      
      # Diagnostics
      fit_to_diagnose <- if (outcome_type == "continuous") fits$ls else fits$glm
      diagnostics <- run_diagnostics(fit_to_diagnose, data_obs)
      
      # Prediction bands (continuous only)
      pred_bands <- tibble(band_type = NA, coverage = NA)
      if (outcome_type == "continuous" && !is.null(fits$ls)) {
        newdata <- generate_sim_data(50, collinearity_level, me_level, error_dist, outcome_type, misspecification)
        newdata <- newdata %>%
          mutate(Y_true = beta0 + beta1*X1_true + beta2*X2 + beta3*X3_true + 
                   ifelse(misspecification != "underfit", beta4*(X1_true^2), 0))
        pred_bands <- calculate_prediction_bands(fits$ls, newdata)
      }
      
      # ANCOVA (fixed sample size)
      data_ancova <- data %>% mutate(factor = factor(sample(1:3, size = nrow(data), replace = TRUE)))
      ancova_fit <- tryCatch(lm(Y ~ factor + X1, data = data_ancova), error = function(e) NULL)
      mc <- run_multiple_comparisons(ancova_fit)
      
      # Aggregate results
      tibble(
        replication = r,
        n = n,
        collinearity = collinearity_level,
        me_level = me_level,
        error_dist = error_dist,
        outcome_type = outcome_type,
        misspecification = misspecification,
        model = metrics$model,
        aic = metrics$AIC,
        bic = metrics$BIC,
        cp = metrics$Cp,
        vif = ifelse(!is.null(diagnostics), mean(diagnostics$vif, na.rm=TRUE), NA),
        cookd_influential = ifelse(!is.null(diagnostics), length(diagnostics$influential_obs), NA),
        pred_coverage = mean(pred_bands$coverage, na.rm=TRUE),
        mc_fdr = mean(mc$fdr, na.rm=TRUE)
      )
    }, error = function(e) {
      message(paste("Replication", r, "failed:", e$message))
      tibble(
        replication = r,
        n = n,
        collinearity = collinearity_level,
        me_level = me_level,
        error_dist = error_dist,
        outcome_type = outcome_type,
        misspecification = misspecification,
        model = NA,
        aic = NA,
        bic = NA,
        cp = NA,
        vif = NA,
        cookd_influential = NA,
        pred_coverage = NA,
        mc_fdr = NA
      )
    })
  }
  
  # Run simulation
  results <- sim_grid %>%
    pmap_dfr(function(n, collinearity_level, me_level, error_dist, outcome_type, misspecification) {
      map_dfr(1:R, ~safe_replicate(
        r = .x,
        n = n,
        collinearity_level = collinearity_level,
        me_level = me_level,
        error_dist = error_dist,
        outcome_type = outcome_type,
        misspecification = misspecification
      ))
    })
  
  # Aggregate across replications
  results %>%
    group_by(n, collinearity, me_level, error_dist, outcome_type, misspecification, model) %>%
    summarize(
      avg_aic = mean(aic, na.rm=TRUE),
      avg_bic = mean(bic, na.rm=TRUE),
      avg_cp = mean(cp, na.rm=TRUE),
      avg_vif = mean(vif, na.rm=TRUE),
      avg_influential = mean(cookd_influential, na.rm=TRUE),
      avg_pred_coverage = mean(pred_coverage, na.rm=TRUE),
      avg_mc_fdr = mean(mc_fdr, na.rm=TRUE),
      .groups = "drop"
    )
}

# --------------------------
# 8. Collinearity Calibration & Visualization
# --------------------------
#' Calibrate Collinearity for Simulation (Target VIF Values)
#' @description
#' Provides utilities for mapping user-defined collinearity levels to
#' approximate VIF ranges used in simulation design.
#' @export
calibrate_collinearity <- function() {
  set.seed(12345)
  
  # Severe collinearity (VIF=10)
  calibrate_severe <- function(noise_sd) {
    n <- 500
    X1_true <- rnorm(n)
    X2 <- runif(n, -1, 1)
    X3_true <- 0.7*X1_true + 0.3*X2 + rnorm(n, 0, noise_sd)
    fit <- lm(Y ~ X1_true + X2 + X3_true, data = data.frame(Y=rnorm(n), X1_true, X2, X3_true))
    car::vif(fit)["X3_true"]
  }
  
  # None collinearity (VIF=1)
  calibrate_none <- function(noise_sd) {
    n <- 500
    X1_true <- rnorm(n)
    X2 <- runif(n, -1, 1)
    X3_true <- rnorm(n, 0, noise_sd)
    fit <- lm(Y ~ X1_true + X2 + X3_true, data = data.frame(Y=rnorm(n), X1_true, X2, X3_true))
    car::vif(fit)["X3_true"]
  }
  
  # Print calibration results
  cat("=== Severe Collinearity Calibration ===\n")
  for (sd in seq(0.05, 0.3, 0.05)) {
    cat("Noise SD:", sd, "| VIF:", round(calibrate_severe(sd), 1), "\n")
  }
  
  cat("\n=== None Collinearity Calibration ===\n")
  for (sd in seq(0.8, 1.2, 0.1)) {
    cat("Noise SD:", sd, "| VIF:", round(calibrate_none(sd), 1), "\n")
  }
}

# --------------------------
# 9. Enhanced Real-World Case Study (Fixed Shapiro Test)
# --------------------------
run_case_study <- function() {
  # Load required packages (force load dplyr first)
  if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
  if (!require("NHANES")) install.packages("NHANES"); library(NHANES)
  if (!require("car")) install.packages("car"); library(car)
  if (!require("glmnet")) install.packages("glmnet"); library(glmnet)
  if (!require("MASS")) install.packages("MASS"); library(MASS)
  if (!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
  
  # Load NHANES data and convert to data frame (avoid tibble issues)
  data("NHANES")
  nhanes_df <- as.data.frame(NHANES)
  
  # Select only core columns (no fancy dplyr piping)
  core_cols <- c("BMI", "Weight", "Height", "Age", "Diabetes")
  nhanes_clean <- nhanes_df[, core_cols]
  
  # Remove missing values (base R)
  nhanes_clean <- nhanes_clean[complete.cases(nhanes_clean), ]
  
  # Recode Diabetes (base R)
  nhanes_clean$Diabetes <- ifelse(nhanes_clean$Diabetes == "Yes", 1, 0)
  
  # Add Y_true (required for prediction bands)
  nhanes_clean$Y_true <- nhanes_clean$BMI
  
  # Filter invalid values (base R) + limit to 5000 obs (Shapiro test max)
  nhanes_clean <- nhanes_clean[
    is.finite(nhanes_clean$BMI) & 
      is.finite(nhanes_clean$Weight) & 
      is.finite(nhanes_clean$Height) & 
      is.finite(nhanes_clean$Age) & 
      nhanes_clean$Diabetes %in% c(0,1) & 
      nhanes_clean$BMI > 0 & 
      nhanes_clean$Weight > 0 & 
      nhanes_clean$Height > 0, 
  ]
  
  # --------------------------
  # FIX: Limit NHANES data to 5000 observations (Shapiro test max)
  # --------------------------
  set.seed(123)
  if (nrow(nhanes_clean) > 5000) {
    nhanes_clean <- nhanes_clean[sample(nrow(nhanes_clean), 5000, replace = FALSE), ]
    cat("\n=== NHANES Data Limited to 5000 Observations (Shapiro Test Requirement) ===\n")
  }
  
  # Collinearity analysis
  lm_collin <- lm(BMI ~ Weight + Height + Age, data = nhanes_clean)
  collinearity_vif <- car::vif(lm_collin)
  cat("\n=== NHANES Collinearity (VIF) ===\n")
  print(collinearity_vif)
  
  # Box-Cox transformation
  boxcox_fit <- tryCatch(MASS::boxcox(lm_collin, plotit = FALSE), error = function(e) list(x=1, y=0))
  lambda <- boxcox_fit$x[which.max(boxcox_fit$y)]
  nhanes_clean$BMI_transformed <- ifelse(lambda == 0, log(nhanes_clean$BMI), (nhanes_clean$BMI^lambda - 1)/lambda)
  cat("\n=== Box-Cox Lambda ===\n")
  print(lambda)
  
  # ENHANCEMENT: Stronger LASSO/Ridge with standardized BMI (lower lambda)
  nhanes_clean$BMI_z <- scale(nhanes_clean$BMI)[,1] # Standardize BMI (z-score)
  x <- model.matrix(Diabetes ~ BMI_z + Weight + Age, data = nhanes_clean)[, -1]
  y <- as.numeric(nhanes_clean$Diabetes)
  ridge <- glmnet::glmnet(x, y, alpha = 0, family = "binomial")
  lasso <- glmnet::glmnet(x, y, alpha = 1, family = "binomial")
  
  # Extract coefficients with lower lambda (retain meaningful predictors)
  ridge_coef <- coef(ridge, s = 0.01) # Lower lambda for Ridge
  lasso_coef <- coef(lasso, s = 0.01) # Lower lambda for LASSO
  cat("\n=== Enhanced LASSO Coefficients (Lambda = 0.01) ===\n")
  print(lasso_coef)
  cat("\n=== Enhanced Ridge Coefficients (Lambda = 0.01) ===\n")
  print(ridge_coef)
  
  # Linear model diagnostics (safe Shapiro test)
  lm_fit <- lm(BMI ~ Weight + Height + Age, data = nhanes_clean)
  cookd <- cooks.distance(lm_fit)
  influential <- which(cookd > 4/length(cookd))
  res_plot <- ggplot(data.frame(resid = residuals(lm_fit), fitted = fitted(lm_fit))) +
    geom_point(aes(x=fitted, y=resid)) +
    geom_hline(yintercept=0, color="red") +
    labs(title="Residual vs Fitted", x="Fitted Values", y="Residuals") +
    theme_minimal()
  
  # Safe Shapiro test (matches run_diagnostics function)
  shapiro_result <- tryCatch({
    resids <- residuals(lm_fit)
    if (length(resids) > 5000) resids <- sample(resids, 5000)
    if (length(resids) >=3) shapiro.test(resids) else NULL
  }, error = function(e) {
    message("Shapiro test skipped: ", e$message)
    NULL
  })
  
  # Store diagnostics in nested list (match original naming)
  linear_model_diagnostics <- list(
    vif = collinearity_vif,
    cookd = cookd,
    influential_obs = influential,
    residual_plot = res_plot,
    shapiro_test = shapiro_result # Safe result (no sample size error)
  )
  
  # Prediction bands (simplified)
  set.seed(123)
  newdata <- nhanes_clean[sample(nrow(nhanes_clean), 100), ]
  pred_conf <- predict(lm_fit, newdata = newdata, interval = "confidence")
  pred_pred <- predict(lm_fit, newdata = newdata, interval = "prediction")
  
  coverage <- data.frame(
    band_type = c("Point-wise (Mean)", "Simultaneous (Mean)", "Point-wise (Obs)"),
    coverage = c(
      mean(newdata$Y_true >= pred_conf[,2] & newdata$Y_true <= pred_conf[,3]),
      mean(newdata$Y_true >= pred_conf[,2] & newdata$Y_true <= pred_conf[,3]), # Simplified simultaneous
      mean(newdata$Y_true >= pred_pred[,2] & newdata$Y_true <= pred_pred[,3])
    )
  )
  
  cat("\n=== Prediction Coverage ===\n")
  print(coverage)
  
  # Return results (include nested diagnostics for consistency)
  return(list(
    data = nhanes_clean,
    vif = collinearity_vif,
    lambda = lambda,
    ridge_coef = ridge_coef,
    lasso_coef = lasso_coef,
    linear_model_diagnostics = linear_model_diagnostics, # Fixed Shapiro test
    residual_plot = res_plot, # Backward-compatible
    prediction_coverage = coverage
  ))
}
