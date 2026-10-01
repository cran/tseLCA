# tests/testthat/test-regression-v1.R
#
# Numerical regression against tseLCA 1.1.1 (see helper-v1-reference.R).
# Guards the 2.0 refactor: estimates, variances, and fit statistics must not
# change except through documented bug fixes.

ref <- readRDS(test_path("fixtures", "v1_reference.rds"))
data_list <- v1_data()

for (nm in names(v1_configs)) {
  test_that(paste("v1.1.1 regression:", nm), {
    skip_on_cran()
    got <- v1_extract(v1_fit(v1_configs[[nm]], data_list))
    expect_equal(got, ref[[nm]], tolerance = 1e-6)
  })
}
