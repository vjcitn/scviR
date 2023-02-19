

test_that("get_ch12_allsce is list of SCE", {
 x = get_ch12_allsce()
 expect_true(inherits(x, "SimpleList"))
 expect_true(inherits(x[[1]], "SingleCellExperiment"))
})

