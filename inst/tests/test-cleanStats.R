
test_that("cleanStats cleans a matrix properly", {
  x <- matrix(1:10, 10, 10, byrow = FALSE)
  x[3, 1] = NA
  x[1, 10] = NaN
  x[10, 10] = NA
  
  A = synlik:::.clean(x)
  
  expect_equal(A[[1]], 3) 
  expect_identical(A[[2]], as.integer(c(1, 3, 10)))
})