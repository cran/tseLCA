#####################################################
### Simulation Script for the tseLCA package:
### tseLCA: Three-Step Estimation for Latent Class Analysis
### -------------------------------------------------
### By: Sam Lee
### E-Mail: samlee@arizona.edu
#####################################################

# Run this script with:
# source(system.file("examples", "tseLCA_sim.R", package = "tseLCA"))

###################################################
### preliminaries
###################################################

# rm(list = ls())
# gc()
# r_opts <- options(
#   prompt = "R> ",
#   continue = "+  ",
#   width = 77,
#   digits = 4,
#   useFancyQuotes = FALSE,
#   warn = 1
# )

# Loading libraries and installing if unavailable

# Install the development version from GitHub
# if (!require("tseLCA")) {
#   if (!require("pak")) {
#     install.packages("pak")
#   }
#   pak::pak("SamLeeBYU/tseLCA")
# }

library(tseLCA)

###################################################
### generate all simulation conditions
###################################################

output.dir <- tempdir() #"tseLCA_output/simulation"
dataset_path <- file.path(output.dir, "sim_datasets.rds")

if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE)
}

if (!file.exists(dataset_path)) {
  message("Generating simulation data (this will take a few minutes)...")

  datasets <- generate_all_conditions(
    n_rep = 500L,
    base_seed = 06262026L,
    sep_levels = c("low", "mid", "high")
  )
  saveRDS(datasets, file = dataset_path)
  message("Completed data generation and saved datasets to: ", dataset_path)
} else {
  message("Loading existing simulation data from: ", dataset_path)
  datasets <- readRDS(dataset_path)
}

###################################################
### obtain measurement models for all conditions
###################################################

# Pre-computing measurement models for each replicate saves time when
# comparing estimators (BCH vs ML, proportional/modal assignment) because
# lca_step1() is the computational bottleneck.
#
# We call multiLCA to obtain two-step vcov estimates.
# While multiLCA cannot accomodate fixed parameter values as input,
# we create covariate-adjusted predictions for the latent class (passed in as starting values) to obtain
# consistent results for multiLCA in the presence of low-separation.
#
# The alternative would be to call multiLCA a bunch of times and take the model
# with the best log-likelihood (which is what we already do for the measurement model).
# This saves some computation time.
#
# Also note that the warning, "Measurement model still failed to converge even after running more iterations. Consider increasing maxIter.measurement and or measurement.tol"
# may trigger in the low separation cases. This is not an issue (this just means that *one* out of the 20 extra measurement models didn't meet the convergence criterion).
# as long as at least one (and preferably the other 19) models converge, we should settle on a well-converged measurement model for step 1.
#
# A poorly converged measurement model can severely bias two-step and three-step estimates.
#
# Structure mirrors datasets: measurement_models[[scenario]][[sep]][[n]][[rep]]

measurement_path <- file.path(output.dir, "measurement_models.rds")

scenarios <- c("covariate", "distal")
sep_levels <- c("low", "mid", "high")
sample_sizes <- c("500", "1000", "2000")

if (file.exists(measurement_path)) {
  cli::cli_alert_info(
    "Loading existing measurement models from: {measurement_path}"
  )
  measurement_models <- readRDS(measurement_path)
} else {
  measurement_models <- list()
}

# Identify which (sc, sep, nn) conditions still need work
n_rep <- length(datasets[[scenarios[1]]][[sep_levels[1]]][[sample_sizes[1]]])

pending <- list()
for (sc in scenarios) {
  for (sep in sep_levels) {
    for (nn in sample_sizes) {
      existing <- measurement_models[[sc]][[sep]][[nn]]
      n_done <- if (is.null(existing)) {
        0L
      } else {
        sum(!vapply(existing, is.null, logical(1L)))
      }
      if (n_done < n_rep) {
        pending[[length(pending) + 1L]] <- list(
          sc = sc,
          sep = sep,
          nn = nn,
          n_done = n_done
        )
      }
    }
  }
}

if (length(pending) == 0L) {
  cli::cli_alert_success("All conditions complete. Nothing to run.")
} else {
  cli::cli_alert_info(
    "Found {length(pending)} condition(s) with incomplete reps. Resuming..."
  )

  total_remaining <- sum(vapply(
    pending,
    function(p) n_rep - p$n_done,
    integer(1L)
  ))

  cli::cli_progress_bar(
    name = "Fitting measurement models",
    total = total_remaining,
    format = paste0(
      "{cli::pb_name} | {cli::pb_bar} {cli::pb_percent} | ",
      "Rep {cli::pb_current}/{cli::pb_total} | ",
      "Elapsed: {cli::pb_elapsed} | ETA: {cli::pb_eta}"
    )
  )

  # Shared objects needed for covariate warm-start
  mGamma.init <- do.call(rbind, bk2018_params$covariate_params)[, -1L]
  p.xz.sim <- function(Z_mat, params) {
    eta_full <- cbind(0, Z_mat %*% params)
    row_max <- apply(eta_full, 1, max)
    exp_eta <- exp(eta_full - row_max)
    exp_eta / rowSums(exp_eta)
  }

  for (cond in pending) {
    sc <- cond$sc
    sep <- cond$sep
    nn <- cond$nn

    reps_data <- datasets[[sc]][[sep]][[nn]]

    reps_fit <- measurement_models[[sc]][[sep]][[nn]]
    if (is.null(reps_fit)) {
      reps_fit <- vector("list", n_rep)
    }

    for (r in seq_len(n_rep)) {
      if (!is.null(reps_fit[[r]])) {
        next
      }

      #Set a unique seed for reproducibility
      set.seed(
        which(scenarios == sc) *
          1e6 +
          which(sep_levels == sep) * 1e4 +
          as.integer(nn) +
          r
      )

      cli::cli_progress_update(
        status = sprintf(
          "scenario=%-10s sep=%-4s n=%-5s rep=%d/%d",
          sc,
          sep,
          nn,
          r,
          n_rep
        )
      )

      reps_fit[[r]] <- tryCatch(
        {
          m.r <- three_step(
            data = reps_data[[r]],
            Y.names = paste0("Y", 1:6),
            n_classes = 3L,
            maxIter.measurement = 5000,
            iter.measurement = 20L,
            R2.threshold = 0.6,
            verbose = FALSE
          )$measurement_model

          #Two-step starting values are currently not available for distal outcomes
          # so we only compare tseLCA to multiLCA's two-step estimation approach
          if (sc == "covariate") {
            c.fitZ <- fitZ_from_fit0(
              fit0 = m.r$fit0,
              data = reps_data[[r]],
              Y.names = paste0("Y", 1:6),
              Zp.names = "Zp",
              maxIter = 500,
              starting_val = mGamma.init
            )

            Y_mat <- as.matrix(reps_data[[r]][, paste0("Y", 1:6)])
            mPhi.init <- m.r$fit0$mPhi

            pi_adj <- p.xz.sim(cbind(1, reps_data[[r]]$Zp), c.fitZ$mGamma)
            log_lik_items <- Y_mat %*%
              log(mPhi.init) +
              (1 - Y_mat) %*% log(1 - mPhi.init)
            log_W <- log(pi_adj) + log_lik_items
            log_W <- log_W -
              apply(log_W, 1, function(x) {
                m <- max(x)
                m + log(sum(exp(x - m)))
              })
            W_init <- exp(log_W)

            reps_data[[r]]$startval <- apply(W_init, 1, which.max)
            has_all_classes <- length(unique(reps_data[[r]]$startval)) == 3L

            c.r <- multilevLCA::multiLCA(
              data = reps_data[[r]],
              Y = paste0("Y", 1:6),
              iT = 3L,
              Z = "Zp",
              startval = if (has_all_classes) "startval" else NULL,
              extout = TRUE,
              verbose = FALSE
            )

            m.r$fitZ <- c.r
            #So we can check which measurement models failed to converged (if any) after the fact
            m.r$fitZ_converged <- abs(diff(tail(c.r$LLKSeries, 2))) < 1e-8
            m.r$fitZ_iters <- c.r$iter
          }

          m.r
        },
        error = function(e) {
          cli::cli_alert_warning(sprintf(
            "three_step failed: scenario=%s sep=%s n=%s rep=%d: %s",
            sc,
            sep,
            nn,
            r,
            conditionMessage(e)
          ))
          NULL
        }
      )
    }

    measurement_models[[sc]][[sep]][[nn]] <- reps_fit
    saveRDS(measurement_models, file = measurement_path)
    cli::cli_alert_success(
      "Saved: scenario={sc} sep={sep} n={nn} ({n_rep} reps)"
    )
  }

  cli::cli_progress_done()
  cli::cli_alert_success(
    "All conditions complete. Final save: {measurement_path}"
  )
}

###################################################
### sim.cond: evaluate one simulation condition
###################################################
#
# cond: character(3), e.g. c("covariate", "low", "500")
#
# Covariate scenario: 5 estimators
#   modal.ml, modal.bch, prop.ml, prop.bch, two_step
#   Target parameter: Zp:C3 slope (true = 1), index 4 in coef vector
#
# Distal scenario: 4 estimators
#   modal.ml, modal.bch, prop.ml, prop.bch
#   Target parameter: mu_C3 (true = 0), index 3 in coef vector
#
# Returns a data.frame with one row per estimator and columns:
#   estimator   : A string to indicator which estimator was evaluated
#                - modal.ml  : Vermunt (2010) ML correction with modal assignment in step 2
#                - modal.bch : BCH correction with modal assignment in step 2
#                - prop.ml   : Vermunt (2010) ML correction with proportional assignment in step 2
#                - prop.bch  : BCH correction with proportional assignment in step 2
#                - two_step  : Bakk & Kuha (2018) two-step estimator (for the covariate scenario only: Zp -> X -> Y)
#   bias        : Monte Carlo mean estimate of the difference between the true value (testing the slope on C2 for covariate estimation and the mu parameter on C3 for the distal outcomes scenario) and the estimated parameter
#   rmse        : Monte Carlo mean estimate of the squared difference between the true values and estimated parameters
#   coverage    : Monte Carlo mean estimate of the coverage of the true parameters (e.g., the proportion that the true values are within estimate +/- 1.96*SE(estimate)) for an alpha-level Wald test
#   se_sd_ratio : Monte Carlo mean estimate of the estimated SE of the estimator divided by the standard deviation of the respective sampled distribution of the corresponding parameter estimates (the closer to 1, the better)
#   n_ok        : How many replications (out of the total 500) for that estimator resulted in a non-degenerate case (happens most frequently for the BCH methods, where the method yields a non-PSD Hessian)

sim.cond <- function(
  datasets,
  measurement_models,
  cond = c("covariate", "low", "500"),
  alpha = 0.05
) {
  sets <- datasets[[cond[1]]][[cond[2]]][[cond[3]]]
  m.mods <- measurement_models[[cond[1]]][[cond[2]]][[cond[3]]]
  R <- length(sets)

  cli::cli_h1(sprintf(
    "Condition: scenario={.val %s}  sep={.val %s}  n={.val %s}  ({R} reps)",
    cond[1],
    cond[2],
    cond[3]
  ))

  # ---- Covariate scenario --------------------------------------------------------------------------------------------------------
  if (cond[1] == "covariate") {
    true_val <- 1 # Zp:C3 slope
    param_idx <- 4L # [Int:C2, Zp:C2, Int:C3, Zp:C3]

    estimators <- c("modal.ml", "modal.bch", "prop.ml", "prop.bch", "two_step")
    n_ok <- matrix(
      FALSE,
      nrow = R,
      ncol = length(estimators),
      dimnames = list(NULL, estimators)
    )
    ests <- matrix(
      NA_real_,
      nrow = R,
      ncol = length(estimators),
      dimnames = list(NULL, estimators)
    )
    ses <- matrix(
      NA_real_,
      nrow = R,
      ncol = length(estimators),
      dimnames = list(NULL, estimators)
    )

    pb <- cli::cli_progress_bar(
      name = "Replicates",
      total = R,
      format = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} | {cli::pb_eta_str} remaining | failures: {.val {sum(!n_ok[seq_len(max(1L, cli::pb_current)), 'modal.ml'])}}"
    )

    for (s in seq_len(R)) {
      dat.s <- sets[[s]]
      m.s <- m.mods[[s]]

      if (is.null(m.s) || isFALSE(m.s$fitZ_converged)) {
        cli::cli_progress_update()
        next
      }

      # ---- modal ML ----------------------------------------------------------------------------------------------------------------
      fit <- tryCatch(
        three_step(
          data = dat.s,
          Y.names = paste0("Y", 1:6),
          Zp.names = "Zp",
          n_classes = 3,
          step1 = m.s,
          use.modal.assignment = TRUE,
          use.bch = FALSE
        ),
        error = function(e) {
          cli::cli_alert_warning("rep {s} modal.ml: {conditionMessage(e)}")
          NULL
        }
      )
      if (!is.null(fit) && !anyNA(fit$three_step)) {
        n_ok[s, "modal.ml"] <- TRUE
        ests[s, "modal.ml"] <- as.vector(fit$three_step)[param_idx]
        ses[s, "modal.ml"] <- sqrt(diag(fit$three_step_vcov))[param_idx]
        if (!is.null(fit$two_step) && !anyNA(fit$two_step)) {
          n_ok[s, "two_step"] <- TRUE
          ests[s, "two_step"] <- as.vector(fit$two_step)[param_idx]
          if (!is.null(fit$two_step_vcov)) {
            ses[s, "two_step"] <- sqrt(diag(fit$two_step_vcov))[param_idx]
          }
        }
      }

      # ---- modal BCH --------------------------------------------------------------------------------------------------------------
      fit <- tryCatch(
        {
          if (cond[2] != "low") {
            three_step(
              data = dat.s,
              Y.names = paste0("Y", 1:6),
              Zp.names = "Zp",
              n_classes = 3,
              step1 = m.s,
              use.bch = TRUE,
              em.maxIter = 500L
            )
          } else {
            NULL
          }
        },
        error = function(e) {
          # cli::cli_alert_warning(sprintf("rep %d: %s", s, conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(fit) && !anyNA(fit$three_step)) {
        n_ok[s, "modal.bch"] <- TRUE
        ests[s, "modal.bch"] <- as.vector(fit$three_step)[param_idx]
        ses[s, "modal.bch"] <- sqrt(diag(fit$three_step_vcov))[param_idx]
      }

      # ---- proportional ML ----------------------------------------------------------------------------------------------------
      fit <- tryCatch(
        three_step(
          data = dat.s,
          Y.names = paste0("Y", 1:6),
          Zp.names = "Zp",
          n_classes = 3,
          step1 = m.s,
          use.modal.assignment = FALSE,
          use.bch = FALSE
        ),
        error = function(e) {
          cli::cli_alert_warning("rep {s} prop.ml: {conditionMessage(e)}")
          NULL
        }
      )
      if (!is.null(fit) && !anyNA(fit$three_step)) {
        n_ok[s, "prop.ml"] <- TRUE
        ests[s, "prop.ml"] <- as.vector(fit$three_step)[param_idx]
        ses[s, "prop.ml"] <- sqrt(diag(fit$three_step_vcov))[param_idx]
      }

      # ---- proportional BCH --------------------------------------------------------------------------------------------------
      fit <- tryCatch(
        {
          if (cond[2] != "low") {
            three_step(
              data = dat.s,
              Y.names = paste0("Y", 1:6),
              Zp.names = "Zp",
              n_classes = 3,
              step1 = m.s,
              use.modal.assignment = FALSE,
              use.bch = TRUE,
              em.maxIter = 500L
            )
          } else {
            NULL
          }
        },
        error = function(e) {
          # cli::cli_alert_warning(sprintf("rep %d: %s", s, conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(fit) && !anyNA(fit$three_step)) {
        n_ok[s, "prop.bch"] <- TRUE
        ests[s, "prop.bch"] <- as.vector(fit$three_step)[param_idx]
        ses[s, "prop.bch"] <- sqrt(diag(fit$three_step_vcov))[param_idx]
      }

      cli::cli_progress_update()
    }

    cli::cli_progress_done()

    # ---- Distal scenario --------------------------------------------------------------------------------------------------------------
  } else if (cond[1] == "distal") {
    true_val <- 0 # mu_C3
    param_idx <- 3L # [mu_C1, mu_C2, mu_C3]

    estimators <- c("modal.ml", "modal.bch", "prop.ml", "prop.bch")
    n_ok <- matrix(
      FALSE,
      nrow = R,
      ncol = length(estimators),
      dimnames = list(NULL, estimators)
    )
    ests <- matrix(
      NA_real_,
      nrow = R,
      ncol = length(estimators),
      dimnames = list(NULL, estimators)
    )
    ses <- matrix(
      NA_real_,
      nrow = R,
      ncol = length(estimators),
      dimnames = list(NULL, estimators)
    )

    pb <- cli::cli_progress_bar(
      name = "Replicates",
      total = R,
      format = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} | {cli::pb_eta_str} remaining"
    )

    for (s in seq_len(R)) {
      dat.s <- sets[[s]]
      m.s <- m.mods[[s]]

      if (is.null(m.s)) {
        cli::cli_progress_update()
        next
      }

      # ---- modal ML ----------------------------------------------------------------------------------------------------------------
      fit <- tryCatch(
        three_step(
          data = dat.s,
          Y.names = paste0("Y", 1:6),
          Zo.name = "Zo",
          n_classes = 3,
          step1 = m.s,
          family = "gaussian",
          use.modal.assignment = TRUE,
          use.bch = FALSE
        ),
        error = function(e) {
          cli::cli_alert_warning("rep {s} modal.ml: {conditionMessage(e)}")
          NULL
        }
      )
      if (!is.null(fit) && !anyNA(fit$three_step)) {
        n_ok[s, "modal.ml"] <- TRUE
        ests[s, "modal.ml"] <- fit$three_step[param_idx]
        ses[s, "modal.ml"] <- sqrt(diag(fit$three_step_vcov))[param_idx]
      }

      # ---- modal BCH --------------------------------------------------------------------------------------------------------------
      fit <- tryCatch(
        {
          if (cond[2] != "low") {
            three_step(
              data = dat.s,
              Y.names = paste0("Y", 1:6),
              Zo.name = "Zo",
              n_classes = 3,
              step1 = m.s,
              use.bch = TRUE,
              em.maxIter = 500L
            )
          } else {
            NULL
          }
        },
        error = function(e) {
          #cli::cli_alert_warning("rep {s} modal.bch: {conditionMessage(e)}")
          NULL
        }
      )
      if (!is.null(fit) && !anyNA(fit$three_step)) {
        n_ok[s, "modal.bch"] <- TRUE
        ests[s, "modal.bch"] <- fit$three_step[param_idx]
        ses[s, "modal.bch"] <- sqrt(diag(fit$three_step_vcov))[param_idx]
      }

      # ---- proportional ML ----------------------------------------------------------------------------------------------------
      fit <- tryCatch(
        three_step(
          data = dat.s,
          Y.names = paste0("Y", 1:6),
          Zo.name = "Zo",
          n_classes = 3,
          step1 = m.s,
          family = "gaussian",
          use.modal.assignment = FALSE,
          use.bch = FALSE
        ),
        error = function(e) {
          cli::cli_alert_warning("rep {s} prop.ml: {conditionMessage(e)}")
          NULL
        }
      )
      if (!is.null(fit) && !anyNA(fit$three_step)) {
        n_ok[s, "prop.ml"] <- TRUE
        ests[s, "prop.ml"] <- fit$three_step[param_idx]
        ses[s, "prop.ml"] <- sqrt(diag(fit$three_step_vcov))[param_idx]
      }

      # ---- proportional BCH --------------------------------------------------------------------------------------------------
      fit <- tryCatch(
        {
          if (cond[2] != "low") {
            three_step(
              data = dat.s,
              Y.names = paste0("Y", 1:6),
              Zo.name = "Zo",
              n_classes = 3,
              step1 = m.s,
              use.modal.assignment = FALSE,
              use.bch = TRUE,
              em.maxIter = 500L
            )
          } else {
            NULL
          }
        },
        error = function(e) {
          #cli::cli_alert_warning("rep {s} prop.bch: {conditionMessage(e)}")
          NULL
        }
      )
      if (!is.null(fit) && !anyNA(fit$three_step)) {
        n_ok[s, "prop.bch"] <- TRUE
        ests[s, "prop.bch"] <- fit$three_step[param_idx]
        ses[s, "prop.bch"] <- sqrt(diag(fit$three_step_vcov))[param_idx]
      }

      cli::cli_progress_update()
    }

    cli::cli_progress_done()
  } else {
    stop("cond[1] must be 'covariate' or 'distal'.", call. = FALSE)
  }

  # ---- Compute metrics ------------------------------------------------------------------------------------------------------------------
  results <- lapply(estimators, function(est) {
    ok <- n_ok[, est]
    e <- ests[ok, est]
    se <- ses[ok, est]
    n <- sum(ok)

    if (n == 0L) {
      cli::cli_alert_danger("{est}: 0 successful replicates")
      return(data.frame(
        estimator = est,
        bias = NA,
        rmse = NA,
        coverage = NA,
        se_sd_ratio = NA,
        n_ok = 0L
      ))
    }

    err <- e - true_val
    bias <- mean(err)
    rmse <- sqrt(mean(err^2))
    coverage <- mean(abs(err) <= qnorm(1 - alpha / 2) * se, na.rm = TRUE)
    se_sd <- mean(se, na.rm = TRUE) / sd(e)

    cli::cli_alert_success(
      "{est}: n_ok={n}  bias={round(bias,4)}  rmse={round(rmse,4)}  cov={round(coverage,3)}  se/sd={round(se_sd,3)}"
    )

    data.frame(
      estimator = est,
      bias = bias,
      rmse = rmse,
      coverage = coverage,
      se_sd_ratio = se_sd,
      n_ok = n
    )
  })

  out <- do.call(rbind, results)
  out$scenario <- cond[1]
  out$separation <- cond[2]
  out$n <- cond[3]
  rownames(out) <- NULL
  out
}

###################################################
### Run all simulation conditions in parallel
###################################################

run_simulation <- function(
  datasets,
  measurement_models,
  n_cores = NULL,
  out_path = NULL
) {
  # ---- Parse available conditions from the datasets structure ----------------
  # Structure: datasets[[scenario]][[separation]][[n]][[rep]]
  conditions <- do.call(
    rbind,
    lapply(
      names(measurement_models),
      function(scenario) {
        do.call(
          rbind,
          lapply(
            names(measurement_models[[scenario]]),
            function(sep) {
              do.call(
                rbind,
                lapply(
                  names(measurement_models[[scenario]][[sep]]),
                  function(n) {
                    data.frame(
                      scenario = scenario,
                      separation = sep,
                      n = n,
                      stringsAsFactors = FALSE
                    )
                  }
                )
              )
            }
          )
        )
      }
    )
  )
  rownames(conditions) <- NULL

  if (is.null(n_cores)) {
    n_cores <- max(c(
      1L,
      min(3 * ((parallel::detectCores() - 1L) %/% 3), nrow(conditions))
    ))
  }

  cli::cli_h1("tseLCA Simulation Study")
  cli::cli_alert_info(sprintf(
    "%d condition(s) detected across %d scenario(s), %d separation level(s), %d sample size(s)",
    nrow(conditions),
    length(unique(conditions$scenario)),
    length(unique(conditions$separation)),
    length(unique(conditions$n))
  ))
  cli::cli_alert_info(sprintf("Parallelising over %d core(s)", n_cores))

  # ---- Run conditions in parallel --------------------------------------------
  if (n_cores > 1L && requireNamespace("parallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(
      cl,
      varlist = c("datasets", "measurement_models", "sim.cond"),
      envir = environment()
    )
    parallel::clusterEvalQ(cl, {
      library(tseLCA)
      library(cli)
    })

    results_list <- parallel::parLapplyLB(
      cl,
      seq_len(nrow(conditions)),
      fun = function(i) {
        cond <- unlist(conditions[i, ])
        tryCatch(
          sim.cond(datasets, measurement_models, cond = cond),
          error = function(e) {
            cli::cli_alert_danger(
              "Condition {cond[1]}/{cond[2]}/{cond[3]} failed: {conditionMessage(e)}"
            )
            data.frame(
              estimator = NA_character_,
              bias = NA_real_,
              rmse = NA_real_,
              coverage = NA_real_,
              se_sd_ratio = NA_real_,
              n_ok = 0L,
              scenario = cond[1],
              separation = cond[2],
              n = cond[3]
            )
          }
        )
      }
    )
  } else {
    results_list <- lapply(
      seq_len(nrow(conditions)),
      function(i) {
        cond <- unlist(conditions[i, ])
        tryCatch(
          sim.cond(datasets, measurement_models, cond = cond),
          error = function(e) {
            cli::cli_alert_danger(
              "Condition {cond[1]}/{cond[2]}/{cond[3]} failed: {conditionMessage(e)}"
            )
            data.frame(
              estimator = NA_character_,
              bias = NA_real_,
              rmse = NA_real_,
              coverage = NA_real_,
              se_sd_ratio = NA_real_,
              n_ok = 0L,
              scenario = cond[1],
              separation = cond[2],
              n = cond[3]
            )
          }
        )
      }
    )
  }

  results <- do.call(rbind, results_list)
  rownames(results) <- NULL

  # ---- Save if requested -----------------------------------------------------
  if (!is.null(out_path)) {
    if (!dir.exists(dirname(out_path))) {
      dir.create(dirname(out_path), recursive = TRUE)
    }
    saveRDS(results, file = out_path)
    cli::cli_alert_success("Results saved to: {.path {out_path}}")
  }

  cli::cli_h2("Simulation complete")
  results
}

sim.results <- run_simulation(
  datasets,
  measurement_models,
  out_path = file.path(output.dir, "sim_results.rds"),
  #Just run sequentially
  n_cores = 1
)

print(sim.results)
