context("internal")

## test_that("range distance works as expected", {
##     vecs <- matrix(runif(10*50), nrow = 50, ncol = 10)
##     expect_equal(simdist:::C_dense_range_dist("COSINE", vecs, 1:ncol(vecs) - 1L, 1:ncol(vecs) - 1L), 
##                  1 - t(vecs) %*% vecs)
## })
