# tests/testthat/test-helpers.R
# Math helpers -- no external deps, all run on CRAN.

test_that("expand_Y one-hot encodes binary items correctly", {
  mY <- matrix(c(0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L), ncol = 2L)
  ivItem <- c(2L, 2L)
  out <- tseLCA:::expand_Y(mY, ivItem)

  expect_equal(dim(out), c(4L, 4L))
  expect_true(all(out %in% c(0L, 1L)))
  #Each item block (pair of columns) sums to 1 per row
  expect_true(all(rowSums(out[, 1:2]) == 1L))
  expect_true(all(rowSums(out[, 3:4]) == 1L))
  #Y=0 column is complement of Y=1 column
  expect_equal(out[, 1L], 1L - out[, 2L])
  expect_equal(out[, 3L], 1L - out[, 4L])
})

test_that("expand_Y handles a 3-category polytomous item", {
  mY <- matrix(c(0L, 1L, 2L, 1L), ncol = 1L)
  ivItem <- c(3L)
  out <- tseLCA:::expand_Y(mY, ivItem)

  expect_equal(dim(out), c(4L, 3L))
  expect_true(all(rowSums(out) == 1L))
  expect_equal(out[1L, ], c(1L, 0L, 0L))
  expect_equal(out[2L, ], c(0L, 1L, 0L))
  expect_equal(out[3L, ], c(0L, 0L, 1L))
})

test_that("expand_Phi (three_step internal) expands n_free x T to K x T", {
  # expand_Phi in three_step.R always treats input as n_free x T:
  # for binary items, row h is P(Y=1|class t); it prepends P(Y=0) = 1-P(Y=1).
  # Two binary items, two classes: n_free=2, T=2 => K=4
  phi_free <- matrix(
    c(
      0.8,
      0.7, # P(Y1=1|C1), P(Y1=1|C2)
      0.3,
      0.6
    ), # P(Y2=1|C1), P(Y2=1|C2)
    nrow = 2L,
    byrow = TRUE
  )
  ivItem <- c(2L, 2L)
  out <- tseLCA:::expand_Phi(phi_free, ivItem)

  expect_equal(dim(out), c(4L, 2L))

  expect_equal(out[1L, ], 1 - phi_free[1L, ]) # P(Y1=0)
  expect_equal(out[2L, ], phi_free[1L, ]) # P(Y1=1)
  expect_equal(out[3L, ], 1 - phi_free[2L, ]) # P(Y2=0)
  expect_equal(out[4L, ], phi_free[2L, ]) # P(Y2=1)
  #Each item's two rows sum to 1
  expect_equal(out[1L, ] + out[2L, ], c(1, 1))
  expect_equal(out[3L, ] + out[4L, ], c(1, 1))
})

test_that("log_lik_matrix returns correct dimensions and non-positive values", {
  set.seed(1L)
  N <- 5L
  K <- 4L
  T <- 2L
  Y <- matrix(sample(0:1, N * K, replace = TRUE), N, K)
  mPhi <- matrix(runif(K * T, 0.1, 0.9), K, T)

  out <- tseLCA:::log_lik_matrix(Y, mPhi)
  expect_equal(dim(out), c(N, T))
  expect_true(all(is.finite(out)))
  expect_true(all(out <= 0)) # log probabilities
})

test_that("log_lik_matrix with mDesign reduces log-likelihood", {
  # When we zero out a column via mDesign, that item no longer contributes
  # to the log-likelihood. With 4 columns and mDesign=0 on columns 3-4,
  # the result equals log_lik_matrix on only the first 2 columns.
  N <- 4L
  K <- 4L
  T <- 2L
  set.seed(2L)
  Y <- matrix(sample(0:1, N * K, replace = TRUE), N, K)
  mPhi <- matrix(runif(K * T, 0.2, 0.8), K, T)

  mDes_full <- matrix(1L, N, K)
  mDes_partial <- mDes_full
  #mask last two columns
  mDes_partial[, 3:4] <- 0L

  out_full <- tseLCA:::log_lik_matrix(Y, mPhi, mDes_full)
  out_partial <- tseLCA:::log_lik_matrix(Y, mPhi, mDes_partial)
  out_2col <- tseLCA:::log_lik_matrix(Y[, 1:2], mPhi[1:2, ])

  #Masking columns reduces log-likelihood
  expect_true(all(out_partial >= out_full - 1e-12))
  #Masked result should equal log_lik on just the first 2 columns
  expect_equal(out_partial, out_2col, tolerance = 1e-12)
})

test_that("compute_pwx_adj returns valid p.wx_mat", {
  set.seed(1L)
  # 2 classes, 2 binary items, N=6; one-hot expanded (K=4 columns)
  Y.obs <- matrix(
    c(
      1L,
      0L,
      0L,
      1L, # item 1: Y=1 => col1=1,col2=0; Y=0 => col1=0,col2=1
      0L,
      1L,
      1L,
      0L,
      1L,
      0L,
      0L,
      1L,
      0L,
      1L,
      1L,
      0L,
      1L,
      0L,
      1L,
      0L,
      0L,
      1L,
      0L,
      1L
    ),
    nrow = 6L,
    byrow = TRUE
  )
  fit0 <- list(
    mPhi = matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2L), # n_free x T
    vPi = c(0.5, 0.5)
  )
  ivItemcat <- c(2L, 2L)

  res <- tseLCA:::compute_pwx_adj(
    Y.obs,
    fit0,
    ivItemcat,
    mDesign = NULL,
    use.modal.assignment = TRUE
  )

  #p.wx_mat is T x T and each column sums to 1
  expect_equal(dim(res$p.wx_mat), c(2L, 2L))
  expect_equal(colSums(res$p.wx_mat), c(1, 1), tolerance = 1e-12)
  #w.is rows sum to 1 (hard assignment: one-hot per row)
  expect_equal(dim(res$w.is), c(6L, 2L))
  expect_equal(rowSums(res$w.is), rep(1L, 6L))
  #All values are 0 or 1 for hard assignment
  expect_true(all(res$w.is %in% c(0L, 1L)))
})
