library(data.table)
library(haldensify)
library(sl3)
library(tmle3)
library(tmle3shift)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(caret)

if (!requireNamespace("EValue", quietly = TRUE)) install.packages("EValue")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")

library(EValue)
library(patchwork)

Sys.setenv("OMP_NUM_THREADS" = 1)
Sys.setenv("MKL_NUM_THREADS" = 1)
Sys.setenv("OPENBLAS_NUM_THREADS" = 1)

options(mc.cores = 1)

# =============================================================================
# E-VALUE ESTIMATION FOR BINARY OUTCOME (excess = 0/1)
# =============================================================================
# tmle3shift estimates EY_delta = P(excess = 1 | do(T = T + delta)).
# The Rd estimated below is a risk difference (RD = EY_delta - EY_0).
# For E-values with a binary outcome, use the causal risk ratio:
# cRR_delta = EY_delta / EY_0.

compute_evalues_binary_tmle <- function(main_results,
                                        stratum_label = "Stratum") {
  required_cols <- c(
    "labs", "EY", "Lower", "Upper", "Rd", "Rd_Lower", "Rd_Upper"
  )
  missing_cols <- setdiff(required_cols, names(main_results))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "[%s] main_results is missing columns: %s",
      stratum_label, paste(missing_cols, collapse = ", ")
    ))
  }
  
  ey_vals <- main_results$EY
  if (any(!is.na(ey_vals) & (ey_vals < 0 | ey_vals > 1))) {
    warning(sprintf(
      "[%s] EY values outside [0, 1]. Verify that 'excess' is binary 0/1.",
      stratum_label
    ))
  }
  
  EY_0 <- main_results[1, "EY"]
  lo_0 <- main_results[1, "Lower"]
  hi_0 <- main_results[1, "Upper"]
  
  if (any(is.na(c(EY_0, lo_0, hi_0))) || EY_0 <= 0) {
    stop(sprintf(
      "[%s] Observed EY is NA or <= 0. cRR and E-values cannot be computed.",
      stratum_label
    ))
  }
  
  SE_0 <- (hi_0 - lo_0) / (2 * 1.96)
  
  cat(sprintf(
    "\n[Stratum: %s]\n  P(excess=1 | observed T) = %.4f [%.4f, %.4f]\n",
    stratum_label, EY_0, lo_0, hi_0
  ))
  
  ev_table <- data.frame(
    Stratum = stratum_label,
    Intervention = main_results$labs,
    EY = main_results$EY,
    EY_Lower95 = main_results$Lower,
    EY_Upper95 = main_results$Upper,
    RD = main_results$Rd,
    RD_Lower95 = main_results$Rd_Lower,
    RD_Upper95 = main_results$Rd_Upper,
    cRR = NA_real_,
    cRR_Lower95 = NA_real_,
    cRR_Upper95 = NA_real_,
    SE_log_cRR = NA_real_,
    EValue_Est = NA_real_,
    EValue_CI = NA_real_,
    stringsAsFactors = FALSE
  )
  
  ev_table[1, c("cRR", "cRR_Lower95", "cRR_Upper95")] <- c(1, 1, 1)
  ev_table[1, c("EValue_Est", "EValue_CI")] <- c(1, 1)

  rr_to_evalue <- function(rr) {
    rr_for_evalue <- ifelse(rr < 1, 1 / rr, rr)
    rr_for_evalue + sqrt(rr_for_evalue * (rr_for_evalue - 1))
  }

  ci_bound_to_evalue <- function(lo, hi) {
    if (lo <= 1 && hi >= 1) return(1)
    if (lo > 1) return(rr_to_evalue(lo))
    if (hi < 1) return(rr_to_evalue(hi))
    return(NA_real_)
  }
  
  for (i in 2:nrow(main_results)) {
    EY_i <- main_results[i, "EY"]
    lo_i <- main_results[i, "Lower"]
    hi_i <- main_results[i, "Upper"]
    lab_i <- main_results[i, "labs"]
    
    if (any(is.na(c(EY_i, lo_i, hi_i)))) {
      warning(sprintf(
        "[%s] Row %d (%s): NA values, skipped.",
        stratum_label, i, lab_i
      ))
      next
    }
    
    if (EY_i <= 0) {
      warning(sprintf(
        "[%s] Row %d (%s): shifted EY <= 0, cRR undefined.",
        stratum_label, i, lab_i
      ))
      next
    }
    
    SE_i <- (hi_i - lo_i) / (2 * 1.96)
    
    cRR_i <- EY_i / EY_0
    SE_log_i <- sqrt((SE_i / EY_i)^2 + (SE_0 / EY_0)^2)
    cRR_lo_i <- exp(log(cRR_i) - 1.96 * SE_log_i)
    cRR_hi_i <- exp(log(cRR_i) + 1.96 * SE_log_i)
    
    ev_table[i, "cRR"] <- cRR_i
    ev_table[i, "cRR_Lower95"] <- cRR_lo_i
    ev_table[i, "cRR_Upper95"] <- cRR_hi_i
    ev_table[i, "SE_log_cRR"] <- SE_log_i
    ev_table[i, "EValue_Est"] <- rr_to_evalue(cRR_i)
    ev_table[i, "EValue_CI"] <- ci_bound_to_evalue(cRR_lo_i, cRR_hi_i)
    
    ev_result <- tryCatch({
      EValue::evalues.RR(est = cRR_i, lo = cRR_lo_i, hi = cRR_hi_i, true = 1)
    }, error = function(e) {
      warning(sprintf(
        "[%s | %s] evalues.RR error: %s",
        stratum_label, lab_i, conditionMessage(e)
      ))
      return(NULL)
    })
    
    if (!is.null(ev_result)) {
      cat(sprintf("\n  -- evalues.RR | %s | %s --\n", stratum_label, lab_i))
      print(ev_result)
    }
  }
  
  cat(sprintf("\nRESULTADOS E-VALUE | Estrato: %s\n", stratum_label))
  print(
    ev_table %>%
      select(
        Intervention, EY, RD, RD_Lower95, RD_Upper95,
        cRR, cRR_Lower95, cRR_Upper95, EValue_Est, EValue_CI
      ) %>%
      mutate(across(where(is.numeric), ~ round(.x, 4))),
    row.names = FALSE
  )
  
  return(ev_table)
}

# --- Data Loading ---
file_path <- "D:/clases/UDES/articulo dengue/stocastic_mil_metros/ci/data_final.csv"

data_all_dengue <- read.csv(file_path, fileEncoding = "latin1")

data_all_dengue <- na.omit(data_all_dengue)

data_all_dengue <- data_all_dengue[data_all_dengue$Year >= 2007, ]

columnas_to_drop <- c('Year_month', 'Period', 'expected',
                      'SIR_mean', 'SIR_lwr95', 'SIR_upr95')

#############################################################################

dataset_m1000 <- data_all_dengue[data_all_dengue$altitude <= 1000, ]

# 1) Calculate total cases per municipality
cases_municipio <- dataset_m1000 %>%
  group_by(DANE) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE))

# 2) Select the 100 municipalities with the most cases
municipios_sel <- cases_municipio %>%
  arrange(desc(total_cases)) %>%
  slice_head(n = 100) 

# 3) Filter original dataset with selected municipalities
dataset_m1000_sel100 <- dataset_m1000 %>%
  semi_join(municipios_sel, by = "DANE")

##############################################################################################################

dataset_m1000_sel100 <- dataset_m1000_sel100 %>%
  select(-all_of(columnas_to_drop))

vars_to_scale <- c('SST12','SST3','SST34','SST4','NATL','SATL','TROP',
                   'SOI','ESOI','cpolr','Zonal_Winds','epac','qbo_u30',
                   'MPI','pop_density','rainfall')

dataset_m1000_sel100[vars_to_scale] <- lapply(
  dataset_m1000_sel100[vars_to_scale],
  function(x) {
    as.numeric(x > median(x, na.rm = TRUE))
  }
)

dataset_m1000_sel100$cosin_month <- cos(2 * pi * dataset_m1000_sel100$Month / 12)
dataset_m1000_sel100$sin_month   <- sin(2 * pi * dataset_m1000_sel100$Month / 12)

dataset_m1000_sel100$Year <- dataset_m1000_sel100$Year - 2007

data_std_dengue <- dataset_m1000_sel100

data_std_dengue <- data_std_dengue %>%
  arrange(DANE, DANE_year, DANE_period)

data_std_dengue <- na.omit(data_std_dengue)

dataset <- data_std_dengue %>%
  select(
    DANE, DANE_year, DANE_period,
    SST12, SST3, SST34, SST4, NATL, SATL, TROP, SOI, # ESOI, #cpolr, Zonal_Winds, epac, qbo_u30,
    rainfall,
    episodes, altitude,
    sin_month, cosin_month, Year,
    temperature, excess
  )


# Apply temperature filter
dataset <- dataset[dataset$temperature >= 15 & dataset$temperature <= 30, ]

# learners used for conditional mean of the outcome
bart_lrnr <- Lrnr_dbarts$new(seed = 12345,nthread = 1)
earth_lrnr <- Lrnr_earth$new(seed = 12345,nthread = 1)
rf_lrnr <- Lrnr_ranger$new(seed = 12345,num.threads = 1)
xgb_lrnr <- Lrnr_xgboost$new(nrounds = 100,max_depth = 5,eta = 0.05,seed = 12345,nthread = 1)


# SL for the outcome regression
#folds <- origami::make_folds(dataset, V = 5)
sl_reg_lrnr <- Lrnr_sl$new(
  learners = list(bart_lrnr, earth_lrnr, rf_lrnr, xgb_lrnr),
  metalearner = Lrnr_nnls$new())

sl3_list_learners("density")

# learners used for conditional densities for (g_n)
haldensify_lrnr <- Lrnr_haldensify$new(
  n_bins = c(5, 10, 20),
  lambda_seq = exp(seq(-1, -10, length = 200))
)
# semiparametric density estimator with homoscedastic errors (HOSE)
hose_rf_lrnr <- make_learner(
  Lrnr_density_semiparametric,
  mean_learner = rf_lrnr
)
# semiparametric density estimator with heteroscedastic errors (HESE)
hese_rf_glm_lrnr <- make_learner(Lrnr_density_semiparametric,
                                 mean_learner = rf_lrnr,
                                 var_learner = earth_lrnr
)

# SL for the conditional treatment density
sl_dens_lrnr <- Lrnr_sl$new(
  learners = list(hose_rf_lrnr, hese_rf_glm_lrnr),
  metalearner = Lrnr_solnp_density$new()
)

learner_list <- list(Y = sl_reg_lrnr, A = sl_dens_lrnr)

# --- Node list ---
print(names(dataset))

node_list <- list(
  W = names(dataset)[4:17], # Explicitly use column names
  A = "temperature",        # Must match exactly
  Y = "excess"              # Must match exactly
)

# --- Ensure Node Variables are Present ---
missing_nodes <- setdiff(c(node_list$A, node_list$Y, node_list$W), names(dataset))
if (length(missing_nodes) > 0) {
  stop(paste("The following nodes are missing from dataset:", paste(missing_nodes, collapse = ", ")))
}

# --- Function to fit TMLE and store results ---
fit_tmle_and_store <- function(shift_val, index) {
  cat("  -> Fitting TMLE for shift =", shift_val, "...\n")
  tmle_spec <- tmle_shift(shift_val = shift_val, shift_fxn = shift_additive, shift_fxn_inv = shift_additive_inv)
  
  tmle_fit <- tryCatch({
    tmle3(tmle_spec, dataset, node_list, learner_list)
  }, error = function(e) {
    cat("  -> ERROR fitting shift =", shift_val, ":", conditionMessage(e), "\n")
    # Print problematic columns if error is related to data types
    if (grepl("numeric|matrix", conditionMessage(e))) {
      cat("     -> Suspect data type issue. Checking column classes:\n")
      problematic_cols <- c(node_list$W[1:min(5, length(node_list$W))], node_list$A, node_list$Y)
      for(col in problematic_cols) {
        cat("        ", col, ":", class(dataset[[col]]), "\n")
      }
    }
    return(NULL)
  })
  
  if (!is.null(tmle_fit)) {
    main_results[index, "EY"] <- tmle_fit$summary$tmle_est
    main_results[index, "Lower"] <- tmle_fit$summary$lower
    main_results[index, "Upper"] <- tmle_fit$summary$upper
  } else {
    main_results[index, "EY"] <- NA
    main_results[index, "Lower"] <- NA
    main_results[index, "Upper"] <- NA
  }
  return(tmle_fit)
}

# --- Dataframe for saving results ---
main_results <- data.frame(
  `labs` = c("Observed temperature", "temperature + 0.5 °C", "temperature + 1.0 °C",
             "temperature + 1.5 °C", "temperature + 2.0 °C"),
  EY = numeric(5),
  Lower = numeric(5),
  Upper = numeric(5),
  Rd = numeric(5),
  Rd_Lower = numeric(5),
  Rd_Upper = numeric(5)
)
main_results[1, "Rd"] <- 0
main_results[1, "Rd_Lower"] <- 0
main_results[1, "Rd_Upper"] <- 0

# --- Fit models ---
cat("\n--- Starting TMLE Fits ---\n")
fit_observed <- fit_tmle_and_store(0, 1)
fit_05 <- fit_tmle_and_store(0.5, 2)
fit_10 <- fit_tmle_and_store(1.0, 3)
fit_15 <- fit_tmle_and_store(1.5, 4)
fit_20 <- fit_tmle_and_store(2.0, 5)

# --- Dataframe for saving results ---
main_results <- data.frame(
  `labs` = c("Observed temperature", "temperature + 0.5 °C", "temperature + 1.0 °C",
             "temperature + 1.5 °C", "temperature + 2.0 °C"),
  EY = c(
    if (!is.null(fit_observed)) fit_observed$summary$tmle_est else NA,
    if (!is.null(fit_05)) fit_05$summary$tmle_est else NA,
    if (!is.null(fit_10)) fit_10$summary$tmle_est else NA,
    if (!is.null(fit_15)) fit_15$summary$tmle_est else NA,
    if (!is.null(fit_20)) fit_20$summary$tmle_est else NA
  ),
  Lower = c(
    if (!is.null(fit_observed)) fit_observed$summary$lower else NA,
    if (!is.null(fit_05)) fit_05$summary$lower else NA,
    if (!is.null(fit_10)) fit_10$summary$lower else NA,
    if (!is.null(fit_15)) fit_15$summary$lower else NA,
    if (!is.null(fit_20)) fit_20$summary$lower else NA
  ),
  Upper = c(
    if (!is.null(fit_observed)) fit_observed$summary$upper else NA,
    if (!is.null(fit_05)) fit_05$summary$upper else NA,
    if (!is.null(fit_10)) fit_10$summary$upper else NA,
    if (!is.null(fit_15)) fit_15$summary$upper else NA,
    if (!is.null(fit_20)) fit_20$summary$upper else NA
  ),
  Rd = numeric(5),
  Rd_Lower = numeric(5),
  Rd_Upper = numeric(5)
)

# Rd vs itself is 0
main_results[1, "Rd"] <- 0
main_results[1, "Rd_Lower"] <- 0
main_results[1, "Rd_Upper"] <- 0

# Rd 
calculrd_rd_manual <- function(index_shift, index_obs = 1) {
  if (is.na(main_results[index_shift, "EY"]) || is.na(main_results[index_obs, "EY"])) {
    warning("can not estimate Rd: Fault EY values.")
    return(c(Rd = NA, Lower = NA, Upper = NA))
  }
  
  psi_obs <- main_results[index_obs, "EY"]
  psi_shift <- main_results[index_shift, "EY"]
  
  # Rd as difference
  rd_est <- psi_shift - psi_obs
  
  # Standard errors
  se_obs <- (main_results[index_obs, "Upper"] - main_results[index_obs, "Lower"]) / (2 * 1.96)
  se_shift <- (main_results[index_shift, "Upper"] - main_results[index_shift, "Lower"]) / (2 * 1.96)
  
  var_rd <- se_shift^2 + se_obs^2
  se_rd <- sqrt(var_rd)
  
  # Confidence interval for Rd
  margin_of_error <- 1.96 * se_rd
  rd_ci_lower <- rd_est - margin_of_error
  rd_ci_upper <- rd_est + margin_of_error
  
  return(c(Rd = rd_est, Lower = rd_ci_lower, Upper = rd_ci_upper))
}

# --- Calculate and Store Rds  ---
manual_rd_05 <- calculrd_rd_manual(2, 1) # +0.5 vs Observed
main_results[2, c("Rd", "Rd_Lower", "Rd_Upper")] <- manual_rd_05

manual_rd_10 <- calculrd_rd_manual(3, 1) # +1.0 vs Observed
main_results[3, c("Rd", "Rd_Lower", "Rd_Upper")] <- manual_rd_10

manual_rd_15 <- calculrd_rd_manual(4, 1) # +1.5 vs Observed
main_results[4, c("Rd", "Rd_Lower", "Rd_Upper")] <- manual_rd_15

manual_rd_20 <- calculrd_rd_manual(5, 1) # +2.0 vs Observed
main_results[5, c("Rd", "Rd_Lower", "Rd_Upper")] <- manual_rd_20

# --- Final results ---
print(main_results)

# Preserve stratum-specific objects before they are overwritten by the next stratum
main_results_m1000 <- main_results
dataset_m1000_tmle <- dataset

evalue_results_m1000 <- compute_evalues_binary_tmle(
  main_results = main_results_m1000,
  stratum_label = "m1000: altitude <= 1000 m"
)


# --- Figure 3: Plot Rds ---
rd_plot_data <- main_results[-1, , drop = FALSE] 
rd_plot_data <- rd_plot_data[!is.na(rd_plot_data$Rd), , drop = FALSE]

if (nrow(rd_plot_data) > 0) {
  rd_plot_data$labs <- factor(rd_plot_data$labs, levels = rd_plot_data$labs)
  
  p_rd <- ggplot(rd_plot_data, aes(x = labs, y = Rd)) +
    geom_point(size = 4, color = "darkblue") +
    geom_errorbar(aes(ymin = Rd_Lower, ymax = Rd_Upper), width = 0.2, color = "darkblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(
      title = "a",
      x = "Interventions",
      y = "Risk difference"
    ) +
    # Define a larger base size for all plot text
    theme_minimal(base_size = 14) + 
    theme(
      # Adjust axis text (category labels)
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      
      # Adjust axis titles
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
      
      # Adjust plot title and subtitle
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 16),
      
      # Optional: Ensure space so large text is not clipped
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # Print
  print(p_rd)
  
} else {
  cat("\n!!!No Rds valid. !!!\n")
}









#################################################################################

dataset_p1000 <- data_all_dengue[(data_all_dengue$altitude > 1000 & data_all_dengue$altitude <= 2300), ]

# 1) Calculate total cases per municipality
cases_municipio <- dataset_p1000 %>%
  group_by(DANE) %>%
  summarise(total_cases = sum(cases, na.rm = TRUE))

# 2) Select the 100 municipalities with the most cases
municipios_sel <- cases_municipio %>%
  arrange(desc(total_cases)) %>%
  slice_head(n = 100)
# 3) Filter original dataset with selected municipalities
dataset_p1000_sel100 <- dataset_p1000 %>%
  semi_join(municipios_sel, by = "DANE")

###########################################
dataset_p1000_sel100 <- dataset_p1000_sel100 %>%
  select(-all_of(columnas_to_drop))

vars_to_scale <- c('SST12','SST3','SST34','SST4','NATL','SATL','TROP',
                   'SOI','ESOI','cpolr','Zonal_Winds','epac','qbo_u30',
                   'MPI','pop_density','rainfall')

dataset_p1000_sel100[vars_to_scale] <- lapply(
  dataset_p1000_sel100[vars_to_scale],
  function(x) {
    as.numeric(x > median(x, na.rm = TRUE))
  }
)

dataset_p1000_sel100$cosin_month <- cos(2 * pi * dataset_p1000_sel100$Month / 12)
dataset_p1000_sel100$sin_month   <- sin(2 * pi * dataset_p1000_sel100$Month / 12)

dataset_p1000_sel100$Year <- dataset_p1000_sel100$Year - 2007

data_std_dengue <- dataset_p1000_sel100

data_std_dengue <- data_std_dengue %>%
  arrange(DANE, DANE_year, DANE_period)

data_std_dengue <- na.omit(data_std_dengue)

dataset <- data_std_dengue %>%
  select(
    DANE, DANE_year, DANE_period,
    SST12, SST3, SST34, SST4, NATL, SATL, TROP, SOI, # ESOI, #cpolr, Zonal_Winds, epac, qbo_u30,
    rainfall,
    episodes, altitude,
    sin_month, cosin_month, Year,
    temperature, excess
  )


# Apply temperature filter
dataset <- dataset[dataset$temperature >= 15 & dataset$temperature <= 30, ]

# learners used for conditional mean of the outcome
bart_lrnr <- Lrnr_dbarts$new(seed = 12345,nthread = 1)
earth_lrnr <- Lrnr_earth$new(seed = 12345,nthread = 1)
rf_lrnr <- Lrnr_ranger$new(seed = 12345,num.threads = 1)
xgb_lrnr <- Lrnr_xgboost$new(nrounds = 100,max_depth = 5,eta = 0.05,seed = 12345,nthread = 1)


# SL for the outcome regression
#folds <- origami::make_folds(dataset, V = 5)
sl_reg_lrnr <- Lrnr_sl$new(
  learners = list(bart_lrnr, earth_lrnr, rf_lrnr, xgb_lrnr),
  metalearner = Lrnr_nnls$new())

sl3_list_learners("density")

# learners used for conditional densities for (g_n)
haldensify_lrnr <- Lrnr_haldensify$new(
  n_bins = c(5, 10, 20),
  lambda_seq = exp(seq(-1, -10, length = 200))
)
# semiparametric density estimator with homoscedastic errors (HOSE)
hose_rf_lrnr <- make_learner(
  Lrnr_density_semiparametric,
  mean_learner = rf_lrnr
)
# semiparametric density estimator with heteroscedastic errors (HESE)
hese_rf_glm_lrnr <- make_learner(Lrnr_density_semiparametric,
                                 mean_learner = rf_lrnr,
                                 var_learner = earth_lrnr
)

# SL for the conditional treatment density
sl_dens_lrnr <- Lrnr_sl$new(
  learners = list(hose_rf_lrnr, hese_rf_glm_lrnr),
  metalearner = Lrnr_solnp_density$new()
)

learner_list <- list(Y = sl_reg_lrnr, A = sl_dens_lrnr)

# --- Node list ---
print(names(dataset))

node_list <- list(
  W = names(dataset)[4:17], # Explicitly use column names
  A = "temperature",        # Must match exactly
  Y = "excess"              # Must match exactly
)

# --- Ensure Node Variables are Present ---
missing_nodes <- setdiff(c(node_list$A, node_list$Y, node_list$W), names(dataset))
if (length(missing_nodes) > 0) {
  stop(paste("The following nodes are missing from dataset:", paste(missing_nodes, collapse = ", ")))
}

# --- Function to fit TMLE and store results ---
fit_tmle_and_store <- function(shift_val, index) {
  cat("  -> Fitting TMLE for shift =", shift_val, "...\n")
  tmle_spec <- tmle_shift(shift_val = shift_val, shift_fxn = shift_additive, shift_fxn_inv = shift_additive_inv)
  
  tmle_fit <- tryCatch({
    tmle3(tmle_spec, dataset, node_list, learner_list)
  }, error = function(e) {
    cat("  -> ERROR fitting shift =", shift_val, ":", conditionMessage(e), "\n")
    # Print problematic columns if error is related to data types
    if (grepl("numeric|matrix", conditionMessage(e))) {
      cat("     -> Suspect data type issue. Checking column classes:\n")
      problematic_cols <- c(node_list$W[1:min(5, length(node_list$W))], node_list$A, node_list$Y)
      for(col in problematic_cols) {
        cat("        ", col, ":", class(dataset[[col]]), "\n")
      }
    }
    return(NULL)
  })
  
  if (!is.null(tmle_fit)) {
    main_results[index, "EY"] <- tmle_fit$summary$tmle_est
    main_results[index, "Lower"] <- tmle_fit$summary$lower
    main_results[index, "Upper"] <- tmle_fit$summary$upper
  } else {
    main_results[index, "EY"] <- NA
    main_results[index, "Lower"] <- NA
    main_results[index, "Upper"] <- NA
  }
  return(tmle_fit)
}

# --- Dataframe for saving results ---
main_results <- data.frame(
  `labs` = c("Observed temperature", "temperature + 0.5 °C", "temperature + 1.0 °C",
             "temperature + 1.5 °C", "temperature + 2.0 °C"),
  EY = numeric(5),
  Lower = numeric(5),
  Upper = numeric(5),
  Rd = numeric(5),
  Rd_Lower = numeric(5),
  Rd_Upper = numeric(5)
)
main_results[1, "Rd"] <- 0
main_results[1, "Rd_Lower"] <- 0
main_results[1, "Rd_Upper"] <- 0

# --- Fit models ---
cat("\n--- Starting TMLE Fits ---\n")
fit_observed <- fit_tmle_and_store(0, 1)
fit_05 <- fit_tmle_and_store(0.5, 2)
fit_10 <- fit_tmle_and_store(1.0, 3)
fit_15 <- fit_tmle_and_store(1.5, 4)
fit_20 <- fit_tmle_and_store(2.0, 5)
















# --- Dataframe for saving results ---
main_results <- data.frame(
  `labs` = c("Observed temperature", "temperature + 0.5 °C", "temperature + 1.0 °C",
             "temperature + 1.5 °C", "temperature + 2.0 °C"),
  EY = c(
    if (!is.null(fit_observed)) fit_observed$summary$tmle_est else NA,
    if (!is.null(fit_05)) fit_05$summary$tmle_est else NA,
    if (!is.null(fit_10)) fit_10$summary$tmle_est else NA,
    if (!is.null(fit_15)) fit_15$summary$tmle_est else NA,
    if (!is.null(fit_20)) fit_20$summary$tmle_est else NA
  ),
  Lower = c(
    if (!is.null(fit_observed)) fit_observed$summary$lower else NA,
    if (!is.null(fit_05)) fit_05$summary$lower else NA,
    if (!is.null(fit_10)) fit_10$summary$lower else NA,
    if (!is.null(fit_15)) fit_15$summary$lower else NA,
    if (!is.null(fit_20)) fit_20$summary$lower else NA
  ),
  Upper = c(
    if (!is.null(fit_observed)) fit_observed$summary$upper else NA,
    if (!is.null(fit_05)) fit_05$summary$upper else NA,
    if (!is.null(fit_10)) fit_10$summary$upper else NA,
    if (!is.null(fit_15)) fit_15$summary$upper else NA,
    if (!is.null(fit_20)) fit_20$summary$upper else NA
  ),
  Rd = numeric(5),
  Rd_Lower = numeric(5),
  Rd_Upper = numeric(5)
)

# Rd vs itself is 0
main_results[1, "Rd"] <- 0
main_results[1, "Rd_Lower"] <- 0
main_results[1, "Rd_Upper"] <- 0

# Rd 
calculrd_rd_manual <- function(index_shift, index_obs = 1) {
  if (is.na(main_results[index_shift, "EY"]) || is.na(main_results[index_obs, "EY"])) {
    warning("can not estimate Rd: Fault EY values.")
    return(c(Rd = NA, Lower = NA, Upper = NA))
  }
  
  psi_obs <- main_results[index_obs, "EY"]
  psi_shift <- main_results[index_shift, "EY"]
  
  # Rd as difference
  rd_est <- psi_shift - psi_obs
  
  # Standard errors
  se_obs <- (main_results[index_obs, "Upper"] - main_results[index_obs, "Lower"]) / (2 * 1.96)
  se_shift <- (main_results[index_shift, "Upper"] - main_results[index_shift, "Lower"]) / (2 * 1.96)
  
  var_rd <- se_shift^2 + se_obs^2
  se_rd <- sqrt(var_rd)
  
  # Confidence interval for Rd
  margin_of_error <- 1.96 * se_rd
  rd_ci_lower <- rd_est - margin_of_error
  rd_ci_upper <- rd_est + margin_of_error
  
  return(c(Rd = rd_est, Lower = rd_ci_lower, Upper = rd_ci_upper))
}

# --- Calculate and Store Rds  ---
manual_rd_05 <- calculrd_rd_manual(2, 1) # +0.5 vs Observed
main_results[2, c("Rd", "Rd_Lower", "Rd_Upper")] <- manual_rd_05

manual_rd_10 <- calculrd_rd_manual(3, 1) # +1.0 vs Observed
main_results[3, c("Rd", "Rd_Lower", "Rd_Upper")] <- manual_rd_10

manual_rd_15 <- calculrd_rd_manual(4, 1) # +1.5 vs Observed
main_results[4, c("Rd", "Rd_Lower", "Rd_Upper")] <- manual_rd_15

manual_rd_20 <- calculrd_rd_manual(5, 1) # +2.0 vs Observed
main_results[5, c("Rd", "Rd_Lower", "Rd_Upper")] <- manual_rd_20

# --- Final results ---
print(main_results)

# Preserve stratum-specific objects before final consolidation
main_results_p1000 <- main_results
dataset_p1000_tmle <- dataset

evalue_results_p1000 <- compute_evalues_binary_tmle(
  main_results = main_results_p1000,
  stratum_label = "p1000: 1000 m < altitude <= 2300 m"
)


# --- Figure 3: Plot Rds ---
rd_plot_data <- main_results[-1, , drop = FALSE] 
rd_plot_data <- rd_plot_data[!is.na(rd_plot_data$Rd), , drop = FALSE]

if (nrow(rd_plot_data) > 0) {
  rd_plot_data$labs <- factor(rd_plot_data$labs, levels = rd_plot_data$labs)
  
  p_rd <- ggplot(rd_plot_data, aes(x = labs, y = Rd)) +
    geom_point(size = 4, color = "darkblue") +
    geom_errorbar(aes(ymin = Rd_Lower, ymax = Rd_Upper), width = 0.2, color = "darkblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(
      title = "a",
      x = "Interventions",
      y = "Risk difference"
    ) +
    # Define a larger base size for all plot text
    theme_minimal(base_size = 14) + 
    theme(
      # Adjust axis text (category labels)
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      
      # Adjust axis titles
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
      
      # Adjust plot title and subtitle
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 16),
      
      # Optional: Ensure space so large text is not clipped
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # Print
  print(p_rd)
  
} else {
  cat("\n!!!No Rds valid. !!!\n")
}




# =============================================================================
# CONSOLIDATED E-VALUE TABLE AND FIGURE
# =============================================================================

evalue_consolidated <- dplyr::bind_rows(
  evalue_results_m1000,
  evalue_results_p1000
)

cat("\n\nTABLA CONSOLIDADA - E-VALUES (2 ESTRATOS)\n")
print(
  evalue_consolidated %>%
    select(
      Stratum, Intervention, EY, RD, cRR, cRR_Lower95, cRR_Upper95,
      EValue_Est, EValue_CI
    ) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))),
  row.names = FALSE
)

write.csv(
  evalue_consolidated,
  "D:/clases/UDES/articulo dengue/stocastic_mil_metros/ci/evalue_binary_results_all_strata.csv",
  row.names = FALSE
)
cat("\nArchivo guardado: evalue_binary_results_all_strata.csv\n")


## Fig 5
plot_data <- evalue_consolidated %>%
  filter(Intervention != "Observed temperature") %>%
  filter(!is.na(EValue_Est)) %>%
  mutate(
    Intervention = factor(
      Intervention,
      levels = main_results_m1000$labs[-1]
    ),
    Stratum = factor(
      Stratum,
      levels = c(
        "m1000: altitude <= 1000 m",
        "p1000: 1000 m < altitude <= 2300 m"
      )
    )
  )

if (nrow(plot_data) > 0) {
  
  colores <- c(
    "m1000: altitude <= 1000 m" = "#2166AC",
    "p1000: 1000 m < altitude <= 2300 m" = "#D6604D"
  )
  stratum_labels <- c(
    "m1000: altitude <= 1000 m" = "altitude <= 1000 m",
    "p1000: 1000 m < altitude <= 2300 m" = "1000 m < altitude <= 2300 m"
  )
  
  p_cRR <- ggplot(
    plot_data,
    aes(x = Intervention, y = cRR, color = Stratum, group = Stratum)
  ) +
    geom_point(size = 3.5) +
    geom_line(linewidth = 0.8) +
    geom_errorbar(
      aes(ymin = cRR_Lower95, ymax = cRR_Upper95),
      width = 0.15,
      linewidth = 0.6
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.7
    ) +
    scale_color_manual(values = colores, labels = stratum_labels) +
    labs(
      title = "a",
      #subtitle = "Causal risk ratio (cRR = EY_delta / EY_0)",
      x = NULL,
      y = "cRR",
      color = "Stratum"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(size = 15, face = "bold"),
      legend.position = "none"
    )
  
  p_ev_est <- ggplot(
    plot_data,
    aes(x = Intervention, y = EValue_Est, color = Stratum, group = Stratum)
  ) +
    geom_point(size = 3.5) +
    geom_line(linewidth = 0.8) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.7
    ) +
    scale_color_manual(values = colores, labels = stratum_labels) +
    labs(
      title = "b",
      #subtitle = "E-value for point estimate",
      x = NULL,
      y = "E-value",
      color = "Stratum"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(size = 15, face = "bold"),
      legend.position = "none"
    )
  
  p_ev_ci <- ggplot(
    plot_data,
    aes(x = Intervention, y = EValue_CI, color = Stratum, group = Stratum)
  ) +
    geom_point(size = 3.5, shape = 17) +
    geom_line(linewidth = 0.8, linetype = "dashed") +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.7
    ) +
    scale_color_manual(values = colores, labels = stratum_labels) +
    labs(
      title = "c",
      #subtitle = "E-value for CI bound closest to the null",
      x = "Temperature intervention",
      y = "E-value (CI)",
      color = "Stratum"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(size = 15, face = "bold"),
      legend.position = "none"
    )
  
  legend_data <- data.frame(
    Stratum = factor(names(stratum_labels), levels = names(stratum_labels)),
    label = unname(stratum_labels),
    x = c(2.1, 4.3),
    y = 1
  )

  make_manual_legend <- function(point_shape = 16) {
    ggplot(legend_data, aes(x = x, y = y, color = Stratum)) +
      annotate("text", x = 0.5, y = 1, label = "Stratum", hjust = 0, size = 3.6) +
      geom_point(size = 3.5, shape = point_shape) +
      geom_text(aes(label = label), hjust = 0, nudge_x = 0.13, size = 3.1) +
      scale_color_manual(values = colores, guide = "none") +
      coord_cartesian(xlim = c(0.3, 6.8), ylim = c(0.75, 1.25), clip = "off") +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))
  }

  legend_points <- make_manual_legend(point_shape = 16)
  legend_triangles <- make_manual_legend(point_shape = 17)

  fig_evalue <- ((p_cRR | p_ev_est) / p_ev_ci / (legend_points | legend_triangles)) +
    patchwork::plot_layout(heights = c(1, 1, 0.16))
  
  print(fig_evalue)

  
}
