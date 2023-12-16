
test_that("detailedPCA", {
  load("testData.RData")
  invisible(capture.output(
    res <- detailedPCA(X)
  ))
  load("storedResult.RData")
  expect_equal(res, Y)
})
