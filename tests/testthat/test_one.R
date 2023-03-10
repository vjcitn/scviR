

test_that("getCh12AllSce is list of SCE", {
 x = getCh12AllSce()
 expect_true(inherits(x, "SimpleList"))
 expect_true(inherits(x[[1]], "SingleCellExperiment"))
})

