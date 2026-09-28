# =============================================================================
# logistic_gcomputation.R
# -----------------------------------------------------------------------------
# Parametric logistic g-computation for the STR+TMLE dengue manuscript.
#
# This script replicates EXACTLY the sample construction and model adjustment
# used by the STR+TMLE analysis (see STR_TMLE.R):
#   - Two altitude strata: low (altitude <= 1000 m) and high
#     (1000 m < altitude <= 2300 m).
#   - Within each stratum, the 100 municipalities with the highest total
#     number of dengue cases (2007-2024) are selected.
#   - Covariates are binarized at the stratum-specific median.
#   - Adjustment set (W): SST12, SST3, SST34, SST4, NATL, SATL, TROP, SOI,
#     rainfall, episodes, altitude, sin_month, cosin_month, Year.
#   - Exposure (A): temperature. Outcome (Y): excess (0/1).
#   - Temperature restricted to [15, 30] degC.
#
# Two models are estimated per stratum:
#   M1 (associational reference): multivariable logistic regression of excess
#      on temperature + W. Reported as odds ratio per +1 degC with
#      cluster-robust (by municipality, DANE) 95% confidence interval.
#   M2 (parametric g-computation): the SAME intervention effect targeted by
#      STR+TMLE (sustained additive shift of temperature by delta degC,
#      delta in {0.5, 1.0, 1.5, 2.0}) is estimated from the fitted logistic
#      model:
#         RD(delta)  = (1/n) sum_i [ plogis(eta_i + bA*delta) - plogis(eta_i) ]
#         cRR(delta) = mean(plogis(eta_i + bA*delta)) / mean(plogis(eta_i))
#      with 95% confidence intervals from a cluster (municipality) bootstrap
#      with B replicates (municipalities resampled with replacement).
#      OR(delta) = exp(bA*delta) is reported as an ILLUSTRATION ONLY
#      (odds ratios are not interchangeable with risk ratios).
#
# Usage:
#   Rscript logistic_gcomputation.R [B] [seed]
#     B     number of bootstrap replicates (default: 1000)
#     seed  random seed for the bootstrap (default: 12345)
#
# Outputs (CSV, written to the same folder as this script):
#   effective_sample_summary.csv   sample sizes per stratum
#   logistic_or_per_1c.csv         M1: OR per +1 degC with cluster-robust CI
#   gcomp_rd_by_stratum_delta.csv  M2: RD_gcomp/cRR_gcomp by stratum and delta
#
# Author: Statistics Agent (Buzz Nest)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Configuration
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

n_boot <- 1000L
boot_seed <- 12345L
if (length(args) >= 1L) {
  parsed <- suppressWarnings(as.integer(args[1L]))
  if (!is.na(parsed) && parsed > 0L) n_boot <- parsed
}
if (length(args) >= 2L) {
  parsed <- suppressWarnings(as.integer(args[2L]))
  if (!is.na(parsed)) boot_seed <- parsed
}

ci_dir    <- "D:/clases/UDES/articulo dengue/stocastic_mil_metros/ci"
data_file <- file.path(ci_dir, "data_final.csv")
out_dir   <- file.path(ci_dir, "g-computation")

delta_grid <- c(0.5, 1.0, 1.5, 2.0)
n_top      <- 100L
z_975      <- qnorm(0.975)

drop_cols <- c("Year_month", "Period", "expected",
               "SIR_mean", "SIR_lwr95", "SIR_upr95")

vars_to_binarize <- c("SST12", "SST3", "SST34", "SST4", "NATL", "SATL",
                      "TROP", "SOI", "ESOI", "cpolr", "Zonal_Winds", "epac",
                      "qbo_u30", "MPI", "pop_density", "rainfall")

model_covariates <- c("SST12", "SST3", "SST34", "SST4", "NATL", "SATL",
                      "TROP", "SOI", "rainfall", "episodes", "altitude",
                      "sin_month", "cosin_month", "Year")

stratum_defs <- list(
  low  = list(label = "low (altitude <= 1000 m)",
              keep  = function(a) a <= 1000),
  high = list(label = "high (1000 m < altitude <= 2300 m)",
              keep  = function(a) a > 1000 & a <= 2300)
)

# -----------------------------------------------------------------------------
# 1. Data preparation (replicates STR_TMLE.R exactly)
# -----------------------------------------------------------------------------
cat("Reading data...\n")
data_all <- read.csv(data_file, fileEncoding = "latin1", stringsAsFactors = FALSE)

# Global complete-case restriction, then calendar window (same order as TMLE).
data_all <- na.omit(data_all)
data_all <- data_all[data_all$Year >= 2007, , drop = FALSE]

# Build one analytic dataset per stratum.
build_stratum_dataset <- function(stratum_name, stratum_def) {
  keep_idx <- stratum_def$keep(data_all$altitude)
  ds <- data_all[keep_idx, , drop = FALSE]

  # Total cases per municipality, ties broken by order of first appearance
  # (identical to the stable arrange(desc()) + slice_head() in STR_TMLE.R).
  uniq_dane <- unique(ds$DANE)
  tot_cases <- vapply(split(ds$cases, factor(ds$DANE, levels = uniq_dane)),
                      sum, numeric(1))
  tot_sorted <- sort(tot_cases, decreasing = TRUE, method = "radix")
  n_avail <- length(tot_sorted)
  n_take <- min(n_top, n_avail)
  top_dane <- names(tot_sorted)[seq_len(n_take)]

  ds <- ds[ds$DANE %in% top_dane, , drop = FALSE]

  # Drop columns that are not used downstream.
  ds <- ds[, setdiff(names(ds), drop_cols), drop = FALSE]

  # Binarize climate/geographic indices at the stratum-specific median
  # (median computed on the selected top-100 sample, before any other filter).
  ds[vars_to_binarize] <- lapply(ds[vars_to_binarize], function(x) {
    as.numeric(x > median(x, na.rm = TRUE))
  })

  # Seasonal Fourier terms and re-centered calendar year.
  ds$sin_month   <- sin(2 * pi * ds$Month / 12)
  ds$cosin_month <- cos(2 * pi * ds$Month / 12)
  ds$Year        <- ds$Year - 2007

  # Sort by municipality and time (as in STR_TMLE.R).
  ds <- ds[order(ds$DANE, ds$DANE_year, ds$DANE_period), , drop = FALSE]
  ds <- na.omit(ds)

  # Final adjustment set: the same 14 covariates used by STR+TMLE
  # (MPI and pop_density are binarized but are NOT part of the model).
  ds <- ds[, c("DANE", "DANE_year", "DANE_period",
               model_covariates, "temperature", "excess"), drop = FALSE]

  # Thermal window applied by the main analysis.
  ds <- ds[ds$temperature >= 15 & ds$temperature <= 30, , drop = FALSE]
  rownames(ds) <- NULL

  return(ds)
}

# -----------------------------------------------------------------------------
# 2. Logistic fit and cluster-robust (sandwich) variance
# -----------------------------------------------------------------------------
fit_logistic <- function(X, y) {
  fit <- glm.fit(x = X, y = y, family = binomial())
  if (!fit$converged) {
    warning("glm.fit did not converge.")
  }
  return(fit)
}

# Cluster-robust covariance (Stata-style CR multiplier G/(G-1)) using the
# canonical-logit score contribution (y - mu) * x.
cluster_robust_vcov <- function(X, y, cluster) {
  fit <- fit_logistic(X, y)
  mu  <- fit$fitted.values

  # Bread: (X' W X)^{-1}, with W = diag(mu*(1-mu)) for the canonical logit.
  w     <- mu * (1 - mu)
  bread <- solve(crossprod(X, w * X))

  # Meat: sum over clusters of outer products of cluster score sums.
  score_mat <- (y - mu) * X
  S <- rowsum(score_mat, group = cluster)          # G x p matrix
  meat <- crossprod(S)

  G <- nrow(S)
  vcov <- (G / (G - 1)) * (bread %*% meat %*% bread)
  dimnames(vcov) <- list(colnames(X), colnames(X))
  return(list(vcov = vcov, fit = fit))
}

# -----------------------------------------------------------------------------
# 3. Point estimates: M1 (OR per +1 degC) and M2 (g-computation)
# -----------------------------------------------------------------------------
compute_point_estimates <- function(ds) {
  X <- model.matrix(
    as.formula(paste("~ temperature +",
                     paste(model_covariates, collapse = " + "))),
    data = ds
  )
  y <- ds$excess

  cr <- cluster_robust_vcov(X, y, cluster = ds$DANE)
  beta  <- cr$fit$coefficients
  bA    <- unname(beta[["temperature"]])
  se_bA <- unname(sqrt(cr$vcov["temperature", "temperature"]))

  or_1c <- exp(bA)
  or_lo <- exp(bA - z_975 * se_bA)
  or_hi <- exp(bA + z_975 * se_bA)

  # G-computation of the shift intervention: delta is added to the observed
  # temperature of each municipality-month; all other covariates are held at
  # their observed values.
  eta0 <- as.numeric(X %*% beta)
  p0   <- plogis(eta0)
  EY0  <- mean(p0)

  point <- list(
    or_per_1c = or_1c, or_lo95 = or_lo, or_hi95 = or_hi,
    bA = bA, se_bA = se_bA, EY0_model = EY0, EY0_obs = mean(y)
  )
  for (d in delta_grid) {
    p1  <- plogis(eta0 + bA * d)
    RD  <- mean(p1) - mean(p0)
    cRR <- mean(p1) / mean(p0)
    point[[paste0("RD_", d)]]  <- RD
    point[[paste0("cRR_", d)]] <- cRR
    point[[paste0("OR_", d)]]  <- exp(bA * d)
    point[[paste0("EY_", d)]]  <- mean(p1)
  }
  return(list(point = point, X = X, y = y, cr = cr))
}

# -----------------------------------------------------------------------------
# 4. Cluster bootstrap for M2 (municipalities resampled with replacement)
# -----------------------------------------------------------------------------
run_cluster_bootstrap <- function(ds, X, y) {
  G <- length(unique(ds$DANE))
  cl_rows <- split(seq_len(nrow(ds)), ds$DANE)
  cl_ids  <- seq_along(cl_rows)

  boot_RD  <- matrix(NA_real_, nrow = n_boot, ncol = length(delta_grid))
  boot_cRR <- matrix(NA_real_, nrow = n_boot, ncol = length(delta_grid))
  boot_OR  <- matrix(NA_real_, nrow = n_boot, ncol = length(delta_grid))
  colnames(boot_RD)  <- paste0("delta_", delta_grid)
  colnames(boot_cRR) <- paste0("delta_", delta_grid)
  colnames(boot_OR)  <- paste0("delta_", delta_grid)
  n_fail <- 0L

  set.seed(boot_seed)
  for (b in seq_len(n_boot)) {
    idx_cl <- sample(cl_ids, size = G, replace = TRUE)
    rows   <- unlist(cl_rows[idx_cl], use.names = FALSE)

    fit_b <- tryCatch(fit_logistic(X[rows, , drop = FALSE], y[rows]),
                      error = function(e) NULL)
    if (is.null(fit_b) || !fit_b$converged) {
      n_fail <- n_fail + 1L
      next
    }
    eta0_b <- as.numeric(X[rows, , drop = FALSE] %*% fit_b$coefficients)
    p0_b   <- plogis(eta0_b)
    bA_b   <- unname(fit_b$coefficients[["temperature"]])
    for (j in seq_along(delta_grid)) {
      d     <- delta_grid[j]
      p1_b  <- plogis(eta0_b + bA_b * d)
      boot_RD[b, j]  <- mean(p1_b) - mean(p0_b)
      boot_cRR[b, j] <- mean(p1_b) / mean(p0_b)
      boot_OR[b, j]  <- exp(bA_b * d)
    }
    if (b %% 100L == 0L) {
      cat(sprintf("    bootstrap replicate %d / %d\n", b, n_boot))
      flush.console()
    }
  }
  return(list(RD = boot_RD, cRR = boot_cRR, OR = boot_OR, n_fail = n_fail))
}

# -----------------------------------------------------------------------------
# 5. Per-stratum analysis
# -----------------------------------------------------------------------------
analyze_stratum <- function(stratum_name, stratum_def) {
  cat(sprintf("\n=== Stratum: %s ===\n", stratum_def$label))
  ds <- build_stratum_dataset(stratum_name, stratum_def)

  n_mun  <- length(unique(ds$DANE))
  n_rows <- nrow(ds)
  n_ev   <- sum(ds$excess)
  cat(sprintf("  Analytic sample: %d municipality-months, %d municipalities, %d excess events\n",
              n_rows, n_mun, n_ev))

  res <- compute_point_estimates(ds)
  pt  <- res$point

  cat(sprintf("  M1 (associational): OR per +1 degC = %.4f (95%% CI %.4f - %.4f) [cluster-robust]\n",
              pt$or_per_1c, pt$or_lo95, pt$or_hi95))
  cat(sprintf("  Baseline risk E[Y]: observed %.4f | model %.4f\n",
              pt$EY0_obs, pt$EY0_model))

  cat(sprintf("  Running cluster bootstrap: B = %d, seed = %d\n", n_boot, boot_seed))
  t0 <- proc.time()
  bt <- run_cluster_bootstrap(ds, res$X, res$y)
  elapsed <- (proc.time() - t0)[["elapsed"]]
  cat(sprintf("  Bootstrap finished in %.1f seconds (%d non-converged replicates)\n",
              elapsed, bt$n_fail))

  tab <- data.frame(
    stratum   = stratum_def$label,
    delta_c   = delta_grid,
    EY_shift  = vapply(delta_grid, function(d) pt[[paste0("EY_", d)]], numeric(1)),
    EY0       = pt$EY0_model,
    RD_gcomp  = vapply(delta_grid, function(d) pt[[paste0("RD_", d)]], numeric(1)),
    RD_lo95   = NA_real_, RD_hi95 = NA_real_,
    cRR_gcomp = vapply(delta_grid, function(d) pt[[paste0("cRR_", d)]], numeric(1)),
    cRR_lo95  = NA_real_, cRR_hi95 = NA_real_,
    OR_illustrative      = vapply(delta_grid, function(d) pt[[paste0("OR_", d)]], numeric(1)),
    OR_illustrative_lo95 = NA_real_, OR_illustrative_hi95 = NA_real_,
    stringsAsFactors = FALSE
  )
  for (j in seq_along(delta_grid)) {
    rd_q  <- quantile(bt$RD[, j], probs = c(0.025, 0.975), na.rm = TRUE)
    crr_q <- quantile(bt$cRR[, j], probs = c(0.025, 0.975), na.rm = TRUE)
    or_q  <- quantile(bt$OR[, j], probs = c(0.025, 0.975), na.rm = TRUE)
    tab$RD_lo95[j]  <- rd_q[[1]]
    tab$RD_hi95[j]  <- rd_q[[2]]
    tab$cRR_lo95[j] <- crr_q[[1]]
    tab$cRR_hi95[j] <- crr_q[[2]]
    tab$OR_illustrative_lo95[j] <- or_q[[1]]
    tab$OR_illustrative_hi95[j] <- or_q[[2]]
  }

  list(ds = ds, sample = c(n_municipalities = n_mun,
                           n_municipality_months = n_rows,
                           n_excess_events = n_ev),
       or = data.frame(stratum = stratum_def$label,
                       OR_per_1c = pt$or_per_1c,
                       OR_lo95 = pt$or_lo95, OR_hi95 = pt$or_hi95,
                       logOR = log(pt$or_per_1c), SE_logOR = pt$se_bA,
                       EY0_observed = pt$EY0_obs, EY0_model = pt$EY0_model),
       m2 = tab, n_boot_fail = bt$n_fail)
}

# -----------------------------------------------------------------------------
# 6. Run both strata and write outputs
# -----------------------------------------------------------------------------
results <- list()
sample_rows <- list()
or_rows     <- list()
m2_rows     <- list()

for (sname in names(stratum_defs)) {
  r <- analyze_stratum(sname, stratum_defs[[sname]])
  results[[sname]] <- r
  sample_rows[[sname]] <- data.frame(
    stratum = stratum_defs[[sname]]$label,
    n_municipalities_analytic = r$sample[["n_municipalities"]],
    n_municipality_months_analytic = r$sample[["n_municipality_months"]],
    n_excess_events = r$sample[["n_excess_events"]],
    bootstrap_replicates = n_boot, bootstrap_seed = boot_seed,
    n_bootstrap_nonconverged = r$n_boot_fail
  )
  or_rows[[sname]] <- r$or
  m2_rows[[sname]] <- r$m2
}

sample_tab <- do.call(rbind, sample_rows)
or_tab     <- do.call(rbind, or_rows)
m2_tab     <- do.call(rbind, m2_rows)

write.csv(sample_tab, file.path(out_dir, "effective_sample_summary.csv"),
          row.names = FALSE)
write.csv(or_tab, file.path(out_dir, "logistic_or_per_1c.csv"),
          row.names = FALSE)
write.csv(m2_tab, file.path(out_dir, "gcomp_rd_by_stratum_delta.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 7. Console summary
# -----------------------------------------------------------------------------
cat("\n\n============================================================\n")
cat("SUMMARY OF RESULTS\n")
cat("============================================================\n")
cat("\n-- Effective analytic sample --\n")
print(sample_tab, row.names = FALSE)

cat("\n-- M1: associational logistic OR per +1 degC (cluster-robust 95% CI) --\n")
print(or_tab[, c("stratum", "OR_per_1c", "OR_lo95", "OR_hi95", "SE_logOR")],
      row.names = FALSE)

cat("\n-- M2: parametric g-computation of the shift intervention --\n")
print(m2_tab[, c("stratum", "delta_c", "EY0", "EY_shift", "RD_gcomp",
                 "RD_lo95", "RD_hi95", "cRR_gcomp", "cRR_lo95", "cRR_hi95")],
      row.names = FALSE)

cat("\nOutput files written to:", out_dir, "\n")
cat("  - effective_sample_summary.csv\n")
cat("  - logistic_or_per_1c.csv\n")
cat("  - gcomp_rd_by_stratum_delta.csv\n")
cat("\nDone.\n")
