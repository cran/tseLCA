# tseLCA/R/clean.R
#
# Shared data-cleaning utilities used by three_step() and fitZ_from_fit0().
#
# The central function is clean_data() which takes the raw data.frame and
# returns consistently prepared matrices for Y (expanded one-hot), Z
# (covariate design), Zo (distal outcome), and mDesign (FIML mask).
#
# Missing data handling
# ------------------------
# Steps 1 & 2 (measurement model, posteriors):
#   Use as much Y data as possible.
#   - Drop rows where ALL Y items are missing (completely uninformative).
#   - If incomplete = FALSE, also drop rows where ANY Y is missing.
#   - Missingness in Z or Zo is irrelevant at this stage.
#
# Step 3 covariate:
#   Restrict to rows that (a) passed the Y filter AND (b) have complete Z.
#
# Step 3 distal:
#   Restrict to rows that (a) passed the Y filter AND (b) have complete Zo
#   (and, when covariates are also modeled, complete covariates).

#' Prepare and validate data for tseLCA estimation
#'
#' @keywords internal
#'
#' @param data       A data.frame.
#' @param Y.names    Character vector of item column names.
#' @param Zp.names   Character vector of covariate column names, or `NULL`.
#'   Ignored when `Zp.formula` is given.
#' @param Zo.name    Single distal outcome column name, or `NULL`.
#' @param incomplete Logical. If `TRUE`, use FIML for partially-observed Y.
#' @param include.intercept Logical. Include an intercept in the covariate
#'   design built from `Zp.names`.
#' @param verbose    Logical. Print row-drop messages.
#' @param Zp.formula One-sided formula for the covariate design (e.g.
#'   `~ age + factor(region)`), or `NULL` to build one from `Zp.names`.
#' @param Y.levels   Named list of indicator categories, as returned by
#'   `.recode_indicators()`, when `data` already holds 0-based codes; `NULL`
#'   recodes the indicators here.
#'
#' @return A named list with:
#' \describe{
#'   \item{Y.obs}{N_Y x K expanded one-hot indicator matrix for Steps 1 & 2.}
#'   \item{mDesign}{N_Y x K design/mask matrix (NULL when incomplete = FALSE).}
#'   \item{ivItemcat}{Integer vector of category counts per item.}
#'   \item{Y.levels}{Named list of the categories of each item.}
#'   \item{keep_Y}{Integer indices of rows kept for Steps 1 & 2 (into original n).}
#'   \item{Z_mat}{n_Z x (Q+1) covariate design matrix, or NULL.}
#'   \item{Zp.formula, Z_terms, Z_xlevels}{The covariate formula, its terms,
#'     and the factor levels used, or NULL.}
#'   \item{keep_step3_Z_in_Y}{Positions of Z-complete rows within keep_Y.}
#'   \item{Zo_mat}{N_Zo x 1 distal outcome matrix, or NULL.}
#'   \item{keep_step3_Zo_in_Y}{Positions of Zo-complete rows within keep_Y.}
#'   \item{keep_step3_Zo}{Indices of Zo-complete rows (into original n).}
#'   \item{keep_step3_Zo_in_Z}{With covariates, positions of the distal rows
#'     within the covariate rows (distal rows then also need complete
#'     covariates); otherwise NULL.}
#' }
clean_data <- function(
  data,
  Y.names,
  Zp.names = NULL,
  Zo.name = NULL,
  incomplete = FALSE,
  include.intercept = TRUE,
  verbose = FALSE,
  Zp.formula = NULL,
  Y.levels = NULL
) {
  if (is.null(Y.levels)) {
    rec <- .recode_indicators(data, Y.names)
    data <- rec$data
    Y.levels <- rec$levels
  }
  Y_raw <- as.matrix(data[, Y.names, drop = FALSE])
  ivItemcat <- lengths(Y.levels)[Y.names]

  # ---- Y row filter -----------------------------------------------------------
  all_Y_missing <- rowSums(!is.na(Y_raw)) == 0L
  any_Y_missing <- rowSums(is.na(Y_raw)) > 0L

  drop_Y <- all_Y_missing
  if (!incomplete) {
    drop_Y <- drop_Y | any_Y_missing
  }

  keep_Y <- which(!drop_Y)

  if (any(drop_Y) && verbose) {
    message(sprintf(
      "%d row(s) dropped from measurement/classification steps (missing Y).",
      sum(drop_Y)
    ))
  }

  Y.obs <- Y_raw[keep_Y, , drop = FALSE]

  # ---- Expand Y and build mDesign ---------------------------------------------
  if (incomplete) {
    Y.obs_exp <- expand_Y(Y.obs, ivItemcat)
    mDesign <- (!is.na(Y.obs_exp)) * 1L
    Y.obs_exp[is.na(Y.obs_exp)] <- 0L
  } else {
    Y.obs[is.na(Y.obs)] <- 0L
    Y.obs_exp <- expand_Y(Y.obs, ivItemcat)
    mDesign <- NULL
  }

  # ---- Covariate design matrix ------------------------------------------------
  Zp.formula <- .covariate_formula(Zp.formula, Zp.names, include.intercept)
  Z_mat <- Z_terms <- Z_xlevels <- NULL
  keep_step3_Z <- NULL
  if (!is.null(Zp.formula)) {
    Z_vars <- all.vars(Zp.formula)
    absent <- setdiff(Z_vars, names(data))
    if (length(absent) > 0L) {
      stop(
        "Covariate(s) not found in `data`: ",
        paste(absent, collapse = ", "),
        call. = FALSE
      )
    }
    any_Z_missing_full <- !stats::complete.cases(data[, Z_vars, drop = FALSE])

    drop_step3_Z <- drop_Y | any_Z_missing_full
    keep_step3_Z <- which(!drop_step3_Z)
    keep_step3_Z_in_Y <- match(keep_step3_Z, keep_Y)

    if (sum(any_Z_missing_full[keep_Y]) > 0L && verbose) {
      message(sprintf(
        "%d row(s) excluded from covariate step (missing Z).",
        sum(any_Z_missing_full[keep_Y])
      ))
    }

    # Design on the covariate-complete rows; unused factor levels dropped
    # (as lm() does).
    mf <- stats::model.frame(
      Zp.formula,
      data[keep_step3_Z, , drop = FALSE],
      drop.unused.levels = TRUE
    )
    Z_terms <- stats::terms(mf)
    Z_xlevels <- stats::.getXlevels(Z_terms, mf)
    Z_mat <- stats::model.matrix(Z_terms, mf)
    attr(Z_mat, "assign") <- NULL
    attr(Z_mat, "contrasts") <- NULL
    rownames(Z_mat) <- NULL
    if (ncol(Z_mat) == 0L) {
      stop("The covariate formula has no terms.", call. = FALSE)
    }

    # Check for linear dependence
    Z_rank <- qr(Z_mat)$rank
    if (Z_rank < ncol(Z_mat)) {
      stop(
        sprintf(
          paste0(
            "Covariate design matrix is rank-deficient (rank %d, %d columns). ",
            "Check for perfectly collinear predictors or a redundant intercept."
          ),
          Z_rank,
          ncol(Z_mat)
        ),
        call. = FALSE
      )
    }
  } else {
    keep_step3_Z_in_Y <- seq_along(keep_Y)
  }

  # ---- Distal outcome matrix --------------------------------------------------
  keep_step3_Zo_in_Z <- NULL
  if (!is.null(Zo.name)) {
    if (!Zo.name %in% names(data)) {
      stop(sprintf("Distal outcome `%s` not found in `data`.", Zo.name), call. = FALSE)
    }
    Zo_mat_full <- as.matrix(data[, Zo.name, drop = FALSE])
    # With covariates, the distal model's class prior P(X | Zp) needs them too.
    needed <- c(Zo.name, if (!is.null(Zp.formula)) all.vars(Zp.formula))
    any_Zo_missing_full <- !stats::complete.cases(data[, needed, drop = FALSE])

    drop_step3_Zo <- drop_Y | any_Zo_missing_full
    keep_step3_Zo <- which(!drop_step3_Zo)
    Zo_mat <- Zo_mat_full[keep_step3_Zo, , drop = FALSE]
    keep_step3_Zo_in_Y <- match(keep_step3_Zo, keep_Y)
    if (!is.null(keep_step3_Z)) {
      keep_step3_Zo_in_Z <- match(keep_step3_Zo, keep_step3_Z)
    }

    if (sum(any_Zo_missing_full[keep_Y]) > 0L && verbose) {
      message(sprintf(
        "%d row(s) excluded from distal step (missing %s).",
        sum(any_Zo_missing_full[keep_Y]),
        if (length(needed) > 1L) "Zo or covariates" else "Zo"
      ))
    }
  } else {
    Zo_mat <- NULL
    keep_step3_Zo_in_Y <- seq_along(keep_Y)
  }

  list(
    Y.obs = Y.obs_exp,
    mDesign = mDesign,
    ivItemcat = ivItemcat,
    Y.levels = Y.levels,
    keep_Y = keep_Y,
    Z_mat = Z_mat,
    Zp.formula = Zp.formula,
    Z_terms = Z_terms,
    Z_xlevels = Z_xlevels,
    keep_step3_Z_in_Y = keep_step3_Z_in_Y,
    Zo_mat = Zo_mat,
    keep_step3_Zo_in_Y = keep_step3_Zo_in_Y,
    keep_step3_Zo = if (!is.null(Zo.name)) keep_step3_Zo else integer(0L),
    keep_step3_Zo_in_Z = keep_step3_Zo_in_Z
  )
}

#' Covariate formula from a formula or a vector of column names
#' @noRd
.covariate_formula <- function(Zp.formula, Zp.names, include.intercept) {
  if (!is.null(Zp.formula)) {
    if (!inherits(Zp.formula, "formula") || length(Zp.formula) != 2L) {
      stop("The covariate formula must be one-sided, e.g. `~ x1 + x2`.", call. = FALSE)
    }
    return(Zp.formula)
  }
  if (is.null(Zp.names)) {
    return(NULL)
  }
  stats::reformulate(
    sprintf("`%s`", Zp.names),
    intercept = include.intercept,
    env = globalenv()
  )
}

#' Categories of one indicator
#'
#' Factors keep their level order (unused levels dropped); logical, character,
#' and numeric indicators are sorted.
#' @noRd
.indicator_levels <- function(x, name) {
  lev <- if (is.factor(x)) {
    levels(droplevels(x))
  } else {
    sort(unique(x[!is.na(x)]))
  }
  if (length(lev) < 2L) {
    stop(
      sprintf("Indicator `%s` has fewer than two observed categories.", name),
      call. = FALSE
    )
  }
  lev
}

#' Recode indicators to 0-based integer category codes
#'
#' Indicators may be factors, logicals, character, or numeric codes (any
#' coding, e.g. 1..K). Each is recoded to 0, ..., K-1 following its
#' categories, which are derived from the data unless `levels` (a named list,
#' e.g. from a previously fitted measurement model) is supplied.
#'
#' @return list(data = `data` with recoded indicator columns, levels = named
#'   list of the categories of each indicator).
#' @noRd
.recode_indicators <- function(data, Y.names, levels = NULL) {
  absent <- setdiff(Y.names, names(data))
  if (length(absent) > 0L) {
    stop(
      "Indicator(s) not found in `data`: ",
      paste(absent, collapse = ", "),
      call. = FALSE
    )
  }
  out <- stats::setNames(vector("list", length(Y.names)), Y.names)
  for (nm in Y.names) {
    x <- data[[nm]]
    lev <- if (!is.null(levels)) levels[[nm]] else .indicator_levels(x, nm)
    code <- match(as.character(x), as.character(lev)) - 1L
    unknown <- !is.na(x) & is.na(code)
    if (any(unknown)) {
      stop(
        sprintf(
          "Indicator `%s` has values outside its categories (%s): %s.",
          nm,
          paste(lev, collapse = ", "),
          paste(utils::head(unique(as.character(x[unknown])), 5L), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    data[[nm]] <- code
    out[[nm]] <- lev
  }
  list(data = data, levels = out)
}

#' Parse and validate the rebase argument
#'
#' @param rebase Character like "C2" or integer class index.
#' @param iT      Total number of classes.
#' @return Integer class index (1-based) to use as reference.
#' @keywords internal
parse_rebase <- function(rebase, iT) {
  if (is.character(rebase)) {
    if (!grepl("^C[0-9]+$", rebase)) {
      stop(
        sprintf(
          '`rebase` must be "C1", "C2", ... "C%d" or an integer. Got: "%s".',
          iT,
          rebase
        ),
        call. = FALSE
      )
    }
    idx <- as.integer(sub("^C", "", rebase))
  } else {
    idx <- as.integer(rebase)
  }
  if (idx < 1L || idx > iT) {
    stop(
      sprintf(
        "`rebase` must be between 1 and %d. Got: %d.",
        iT,
        idx
      ),
      call. = FALSE
    )
  }
  idx
}

#' Normalize row/column names of a fitZ$mGamma matrix
#'
#' A plain `multiLCA` object uses `rownames` like `"gamma(Intercept|C)"` and
#' `"gamma(Zp|C)"`. This function strips the `gamma(...)` wrapper so names
#' match the clean format used throughout tseLCA (`"Intercept"`, `"Zp"`, etc.)
#' and ensures column names are `"C2"`, `"C3"`, etc.
#'
#' @param fitZ  A fitZ-like list with at least `$mGamma`.
#' @param Zp.names Character vector of covariate column names (used to set
#'   clean rownames when the raw names can't be parsed). If `NULL`, rownames
#'   are stripped from the `gamma(X|C)` pattern only.
#' @param n_classes Integer. Total number of classes (used to derive clean
#'   column names if they are non-standard).
#' @return `fitZ` with normalized `$mGamma` row/col names.
#' @keywords internal
normalize_fitZ_names <- function(fitZ, Zp.names = NULL, n_classes = NULL) {
  if (is.null(fitZ) || is.null(fitZ$mGamma)) {
    return(fitZ)
  }

  mG <- fitZ$mGamma

  # ---- Normalize rownames ----------------------------------------------------
  rn <- rownames(mG)
  if (!is.null(rn)) {
    # Strip "gamma(X|C)" -> "X"
    rn_clean <- sub("^gamma\\((.+)\\|C\\)$", "\\1", rn)
    # Also handle "gamma(X)" without the |C suffix
    rn_clean <- sub("^gamma\\((.+)\\)$", "\\1", rn_clean)
    # multilevLCA labels the intercept "Intercept"; use R's "(Intercept)"
    rn_clean[rn_clean == "Intercept"] <- "(Intercept)"
    rownames(mG) <- rn_clean
  } else if (!is.null(Zp.names)) {
    rownames(mG) <- c("(Intercept)", Zp.names)
  }

  # ---- Normalize colnames ----------------------------------------------------
  cn <- colnames(mG)
  if (!is.null(cn)) {
    # Already clean ("C2", "C3", ...)
    if (!all(grepl("^C[0-9]+$", cn))) {
      # Non-standard: derive from n_classes if available
      if (!is.null(n_classes)) {
        colnames(mG) <- paste0("C", seq_len(ncol(mG)) + 1L)
      }
    }
  } else if (!is.null(n_classes)) {
    colnames(mG) <- paste0("C", seq_len(ncol(mG)) + 1L)
  }

  fitZ$mGamma <- mG
  fitZ
}

#' Permute class columns of a fitZ object to match a new reference class
#'
#' Rebases a `fitZ` object (output of `fitZ_from_fit0` or
#' `fitZ_from_multiLCA`) so that `ref_idx` becomes the reference class.
#' This involves:
#' \enumerate{
#'   \item Rebasing `$mGamma`: reconstructing the full T-column log-ratio
#'     matrix, subtracting the new reference column, and dropping it.
#'   \item Propagating through `$Varmat_cor` with the delta method: the
#'     rebasing transformation is linear (`gamma_new = A * gamma_old`)
#'     so the vcov transforms exactly as `A %*% V %*% t(A)`.
#'   \item Updating all column names.
#' }
#'
#' @param fitZ    Output of `fitZ_from_fit0()` or `fitZ_from_multiLCA()`.
#' @param ref_idx Integer. New reference class (1-based index into the T
#'   classes as currently ordered in `fitZ`).
#' @return `fitZ` with `$mGamma`, `$Varmat_cor`, and names updated.
#' @keywords internal
permute_fitZ_classes <- function(fitZ, ref_idx) {
  if (is.null(fitZ)) {
    return(NULL)
  }

  # Normalize names first so all downstream logic sees clean "Intercept"/"Zp"
  # rownames and "C2"/"C3" colnames regardless of whether fitZ came from
  # fitZ_from_fit0, fitZ_from_multiLCA, or a raw multiLCA call.
  fitZ <- normalize_fitZ_names(fitZ, n_classes = ncol(fitZ$mGamma) + 1L)

  mGamma <- fitZ$mGamma # Q x (T-1): cols = non-ref classes (C2..CT)
  Q <- nrow(mGamma)
  iT <- ncol(mGamma) + 1L # total number of classes

  if (ref_idx == 1L) {
    return(fitZ)
  }

  # ---- Step 1: rebase mGamma --------------------------------------------------
  # Reconstruct Q x T full log-ratio matrix (column 1 = 0, reference)
  gamma_full <- cbind(0, mGamma) # Q x T

  new_ref_col <- gamma_full[, ref_idx, drop = FALSE] # Q x 1
  gamma_rebased <- gamma_full - as.vector(new_ref_col) # Q x T, col ref_idx = 0

  # Drop the new reference column and keep remaining classes in ascending order
  keep_cols <- seq_len(iT)[-ref_idx] # T-1 indices
  gamma_new <- gamma_rebased[, keep_cols, drop = FALSE]
  colnames(gamma_new) <- paste0("C", keep_cols)
  rownames(gamma_new) <- rownames(mGamma)
  fitZ$mGamma <- gamma_new

  # ---- Step 2: propagate Varmat_cor with the delta method ----------------------
  # The rebasing is a linear map on the vec(mGamma) parameter vector.
  # Stacking columns: theta = vec(mGamma), length = Q*(T-1).
  # After rebasing: theta_new = (A_kron_I_Q) * theta_old
  # where A is the (T-1) x (T-1) contrast matrix acting on class columns.
  #
  # A = gamma_rebased[, keep_cols] expressed as a linear function
  # of gamma_full[, -1] (the original non-ref columns).
  #
  # gamma_full[, keep_cols] = gamma_full[, keep_cols]
  #                         - gamma_full[, ref_idx] * 1'
  # i.e. A_col[j] = e_{keep_cols[j]-1} - e_{ref_idx-1}  (in the T-1 basis)
  # where e_k is the k-th standard basis vector of R^{T-1},
  # with the convention that the "C1" column has index 0 (not in basis).
  #
  # For the original columns j = 2..T (indexed 1..T-1 in mGamma):
  #   A[i, j] = I(keep_cols[i] - 1 == j)
  #           - I(ref_idx - 1   == j)

  if (!is.null(fitZ$Varmat_cor) || !is.null(fitZ$raw_fit$Varmat_cor)) {
    V <- if (!is.null(fitZ$Varmat_cor)) {
      fitZ$Varmat_cor
    } else {
      fitZ$raw_fit$Varmat_cor
    } # Q*(T-1) x Q*(T-1)

    # Build (T-1) x (T-1) column-space contrast matrix A
    # Original columns are indexed 1..(T-1) corresponding to classes C2..CT
    old_non_ref <- seq_len(iT - 1L) # 1..(T-1) indexing into mGamma columns
    # keep_cols are class indices (1-based), need to map to mGamma col indices
    keep_mGamma_cols <- keep_cols - 1L # subtract 1 because C1 is not in mGamma
    ref_mGamma_col <- ref_idx - 1L # column of the new ref in old mGamma

    A <- matrix(0, iT - 1L, iT - 1L)
    for (i in seq_len(iT - 1L)) {
      j_direct <- keep_mGamma_cols[i]
      if (j_direct >= 1L && j_direct <= iT - 1L) {
        A[i, j_direct] <- 1
      }
      if (ref_mGamma_col >= 1L && ref_mGamma_col <= iT - 1L) {
        A[i, ref_mGamma_col] <- A[i, ref_mGamma_col] - 1
      }
    }

    # Full transformation: (A kron I_Q)
    A_kron <- kronecker(A, diag(Q)) # Q*(T-1) x Q*(T-1)
    V_new <- A_kron %*% V %*% t(A_kron)

    # Update names
    param_names <- as.vector(outer(
      rownames(mGamma),
      paste0("C", keep_cols),
      paste,
      sep = ":"
    ))
    rownames(V_new) <- param_names
    colnames(V_new) <- param_names

    if (!is.null(fitZ$Varmat_cor)) {
      fitZ$Varmat_cor <- V_new
    }
    if (!is.null(fitZ$raw_fit$Varmat_cor)) {
      fitZ$raw_fit$Varmat_cor <- V_new
    }
  }

  fitZ
}

#' Permute class columns of a fit0 object so that class ref_idx is first
#'
#' Reorders columns of mPhi and vPi so that the desired reference class
#' becomes column 1 before estimation. This ensures the multinomial logit
#' is parameterized with the correct baseline from the start.
#'
#' @param fit0    Raw multilevLCA fit object (has $mPhi and $vPi).
#' @param ref_idx Integer. Class index to move to position 1.
#' @return fit0 with columns permuted.
#' @keywords internal
permute_fit0_classes <- function(fit0, ref_idx) {
  if (ref_idx == 1L) {
    return(fit0)
  }
  iT <- ncol(fit0$mPhi)
  ord <- c(ref_idx, seq_len(iT)[-ref_idx])
  fit0$mPhi <- fit0$mPhi[, ord, drop = FALSE]
  fit0$vPi <- fit0$vPi[ord]
  fit0
}

#' Extract Y.exp, mDesign, posteriors from a multilevLCA mU matrix
#'
#' `fit0$mU` from \pkg{multilevLCA} stores data already in one-hot expanded
#' form: each item k occupies R_k consecutive columns (one per category),
#' followed by T columns of posterior class probabilities.
#'
#' For dichotomous items (R_k=2) the two columns are stored. For polytomous
#' items (R_k>2) all R_k columns are stored. This function first compresses
#' the expanded Y back to integer codes, then re-expands
#' consistently with \code{expand_Y} so downstream functions receive the correct
#' n x sum(R_k) matrix.
#'
#' @param fit0       Raw multilevLCA fit object with \code{$mU}, \code{$mPhi},
#'   \code{$vPi}.
#' @param ivItemcat  Integer vector of category counts per item (length K).
#'   If \code{NULL}, inferred from \code{fit0$mPhi} dimensions.
#'
#' @return A list with:
#' \describe{
#'   \item{Y.exp}{n x sum(R_k) expanded one-hot matrix (NAs replaced with 0).}
#'   \item{mDesign}{n x sum(R_k) design/mask matrix. \code{NULL} if no missing.}
#'   \item{ivItemcat}{Integer vector of category counts per item.}
#'   \item{u_post}{n x T posterior class probability matrix from \code{mU}.}
#' }
#' @keywords internal
extract_Y_from_mU <- function(fit0, ivItemcat = NULL) {
  mU <- fit0$mU
  if (is.null(mU)) {
    stop(
      "fit0$mU is NULL -- multilevLCA must be run with mU stored. Run multiLCA again with etxout=TRUE.",
      call. = FALSE
    )
  }

  iT <- length(fit0$vPi)

  # ---- Infer ivItemcat from column names if not supplied ---------------------
  # mU column structure:
  #   Dichotomous item (K=2) : 1 column  (raw 0/1)
  #   Polytomous item (K>2)  : K columns (one-hot, suffixed .0/.1/.../.K-1)
  # followed by T posterior columns (C1, C2, ...).
  if (is.null(ivItemcat)) {
    cn <- colnames(mU)
    if (is.null(cn)) {
      stop(
        "ivItemcat must be supplied when fit0$mU has no column names.",
        call. = FALSE
      )
    }
    y_names <- cn[seq_len(ncol(mU) - iT)]
    # Polytomous columns end in ".0", ".1", etc.; dichotomous do not
    has_suffix <- grepl("\\.[0-9]+$", y_names)
    item_base <- ifelse(has_suffix, sub("\\.[0-9]+$", "", y_names), y_names)
    # Count columns per unique item (preserving order)
    unique_items <- unique(item_base)
    ivItemcat <- vapply(
      unique_items,
      \(nm) {
        n <- sum(item_base == nm)
        if (n == 1L) 2L else as.integer(n) # 1 col -> dichotomous (K=2)
      },
      integer(1L)
    )
  }

  # Number of Y columns in mU (dichotomous = 1 col, polytomous = K cols)
  n_mU_Y_cols <- sum(ifelse(ivItemcat == 2L, 1L, ivItemcat))
  mY_raw <- mU[, seq_len(n_mU_Y_cols), drop = FALSE]
  u_post <- mU[, (n_mU_Y_cols + 1L):(n_mU_Y_cols + iT), drop = FALSE]
  mode(u_post) <- "double"

  # ---- Compress to N x H integer matrix --------------------------------------
  # Walk items; dichotomous columns are already 0/1 integer codes.
  # Polytomous columns are one-hot blocks -> compress to 0-based integer code.
  H <- length(ivItemcat)
  mY_int <- matrix(NA_integer_, nrow(mY_raw), H)
  col <- 1L

  for (h in seq_len(H)) {
    K_h <- ivItemcat[h]
    if (K_h == 2L) {
      # Single column, already 0/1
      mY_int[, h] <- as.integer(mY_raw[, col])
      col <- col + 1L
    } else {
      # K_h columns per item. multilevLCA codes these two ways: one-hot
      # (1 = observed category, all NA = missing item) on its FIML path, and
      # 0 for the other categories with NA marking the observed category on
      # its listwise path. Decode both.
      block <- mY_raw[, col:(col + K_h - 1L), drop = FALSE]
      mY_int[, h] <- apply(block, 1L, \(row) {
        if (all(is.na(row))) {
          NA_integer_
        } else if (any(row == 1, na.rm = TRUE)) {
          which(row == 1)[1L] - 1L
        } else if (sum(is.na(row)) == 1L) {
          which(is.na(row)) - 1L
        } else {
          NA_integer_
        }
      })
      col <- col + K_h
    }
  }

  # ---- Re-expand to one-hot with mDesign -------------------------------------
  if (anyNA(mY_int)) {
    Y_exp <- expand_Y(mY_int, ivItemcat)
    mDesign <- (!is.na(Y_exp)) * 1L
    Y_exp[is.na(Y_exp)] <- 0L
  } else {
    Y_exp <- expand_Y(mY_int, ivItemcat)
    mDesign <- NULL
  }

  list(
    Y.exp = Y_exp,
    mDesign = mDesign,
    ivItemcat = ivItemcat,
    u_post = u_post
  )
}

#' Prepare the data for three_step()
#'
#' Validates the distal outcome for its family and encodes it (see
#' .encode_distal()), then runs clean_data() on data whose indicators are
#' already recoded to 0-based codes with categories `Y.levels`. Returns
#' clean_data()'s output plus `data` (with the encoded outcome),
#' `zo_levels`, `Y.names`, `Zp.names`, and `Zo.name`.
#' @noRd
.prepare_data <- function(
  data,
  Y.names,
  Zp.names,
  Zo.name,
  family,
  opts,
  Y.levels,
  Zp.formula = NULL
) {
  zo_levels <- NULL
  if (!is.null(Zo.name)) {
    if (!Zo.name %in% names(data)) {
      stop(sprintf("Distal outcome `%s` not found in `data`.", Zo.name), call. = FALSE)
    }
    enc <- .encode_distal(data[[Zo.name]], Zo.name, family)
    data[[Zo.name]] <- enc$z
    zo_levels <- enc$levels
  }

  cd <- clean_data(
    data = data,
    Y.names = Y.names,
    Zp.names = Zp.names,
    Zo.name = Zo.name,
    incomplete = opts$incomplete,
    include.intercept = opts$include.intercept,
    verbose = opts$verbose,
    Zp.formula = Zp.formula,
    Y.levels = Y.levels
  )
  c(
    cd,
    list(
      data = data,
      zo_levels = zo_levels,
      Y.names = Y.names,
      Zp.names = Zp.names,
      Zo.name = Zo.name
    )
  )
}

#' Validate and encode a distal outcome for its family
#'
#' * gaussian: numeric.
#' * poisson: non-negative whole numbers.
#' * binomial: 0/1 numeric, logical, or a two-category factor/character
#'   (its second category is coded 1).
#' * multinomial: any categorical variable with at least two categories,
#'   coded 1..C (`levels` records the categories).
#' @return list(z = encoded outcome, levels = categories or NULL).
#' @noRd
.encode_distal <- function(z, name, family) {
  obs <- z[!is.na(z)]
  fail <- function(what) {
    stop(sprintf("Distal outcome `%s` must be %s for family = \"%s\".", name, what, family),
         call. = FALSE)
  }
  if (family == "multinomial") {
    f <- factor(z)
    if (nlevels(f) < 2L) fail("categorical with at least 2 distinct categories")
    return(list(z = as.integer(f), levels = levels(f)))
  }
  if (family == "binomial") {
    if (is.logical(z)) {
      return(list(z = as.integer(z), levels = NULL))
    }
    if (is.factor(z) || is.character(z)) {
      f <- droplevels(factor(z))
      if (nlevels(f) != 2L) fail("binary (two categories)")
      return(list(z = as.integer(f) - 1L, levels = levels(f)))
    }
    if (!is.numeric(z) || !all(obs %in% c(0, 1))) fail("binary (0/1, logical, or two categories)")
    return(list(z = z, levels = NULL))
  }
  if (!is.numeric(z)) fail("numeric")
  if (family == "poisson" && !all(obs >= 0 & obs == round(obs))) {
    fail("a count (non-negative whole numbers)")
  }
  list(z = z, levels = NULL)
}
