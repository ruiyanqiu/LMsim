# --------------------------
# Demo Script for LMsim (STAT 5430 Project - Ruiyan Qiu)
# Fully Reproducible & Error-Free
# --------------------------

# === STEP 0: Critical Setup (Auto-Run) ===
# Set working directory (UPDATE THIS PATH TO MATCH YOUR LMsim FOLDER!)
setwd("/Users/qry/Desktop/UVA BSN/2025 Fall/STAT 5430 Stat Comp with Python and R/5430_project/LMsim/")

# Install/load required packages (no manual setup needed)
required_pkgs <- c("tidyverse", "glmnet", "car", "agricolae", "purrr", "NHANES", "MASS", "dplyr", "ggplot2")
lapply(required_pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
})

# Load core functions (from R/LM_function.R)
source("R/LM_function.R")

# === STEP 1: Calibrate Collinearity (Target VIF Values) ===
cat("=== Starting Collinearity Calibration ===\n")
calibrate_collinearity()

# === STEP 2: Run Enhanced Simulation (100 Replicates for Stability) ===
sim_grid <- expand.grid(
  n = 200,
  collinearity_level = "none",
  me_level = "none",
  error_dist = "normal",
  outcome_type = "continuous",
  misspecification = "correct",
  stringsAsFactors = FALSE
)
sim_results <- run_full_simulation(sim_grid, R = 100)

# Filter to LS model (only model with valid AIC)
sim_results_ls <- sim_results %>% filter(model == "LS")
cat("\n=== Enhanced Simulation Results (LS Model Only) ===\n")
print(sim_results_ls[, c("model", "avg_aic", "avg_pred_coverage")])

# === STEP 3: Run Enhanced Real-World Case Study (Error-Free) ===
case_study_results <- run_case_study()
cat("\n=== Real-World Case Study: Prediction Coverage ===\n")
print(case_study_results$prediction_coverage)

# Fixed: Print residual plot (no NULL, no Shapiro error)
cat("\n=== Case Study: Residual vs. Fitted Plot ===\n")
print(case_study_results$linear_model_diagnostics$residual_plot)

# Optional: Print Shapiro test results (if available)
if (!is.null(case_study_results$linear_model_diagnostics$shapiro_test)) {
  cat("\n=== Shapiro-Wilk Normality Test (Residuals) ===\n")
  print(case_study_results$linear_model_diagnostics$shapiro_test)
} else {
  cat("\n=== Shapiro-Wilk Test: Skipped (Sample Size Issue) ===\n")
}

# Final confirmation
cat("\n=== All Analyses Completed Successfully! ===")
