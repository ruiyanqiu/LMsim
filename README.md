# STAT 5430 Project: LMsim (Linear Model Simulation)
Author: Ruiyan Qiu  
Working Directory: `/Users/qry/Desktop/UVA BSN/2025 Fall/STAT 5430 Stat Comp with Python and R/5430_project/LMsim/`

## Folder Structure
- `R/LM_function.R`: Core functions for simulation, model fitting, diagnostics, and case study
- `demo.R`: Reproducible script to run the full analysis

## Overview
This project evaluates linear model performance under collinearity, measurement error, misspecification, and non-normality. Key features:
1. Simulate controlled linear model data
2. Fit LS/GLM/Ridge/LASSO/stepwise models
3. Calculate model selection metrics (AIC/BIC/Mallow's Cp)
4. Run diagnostics (VIF/Cook's distance/safe Shapiro-Wilk test)
5. Evaluate prediction intervals and multiple comparisons
6. Real-world case study with NHANES data

## How to Run
1. Install dependencies: `tidyverse, glmnet, car, agricolae, purrr, NHANES, MASS, dplyr, ggplot2`
2. Source the core function file:
   ```r
   source("R/LM_function.R")  # Run this from the LMsim folder
    ```
3. Example execution: using `demo.R` to run `R/LM_function.R`
