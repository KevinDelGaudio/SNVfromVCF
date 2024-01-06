library(testthat)

test_that("Running TypeCounts with an invalid data.frame", {
    df <- data.frame(col1 = rep(1, 10), col2 = rep(2, 10))
    expect_error(TypeCounts(df, plot = FALSE), "Invalid Mutations Table")
})

test_that("Running TypeCounts with an valid data.frame", {
    data("vcf_chr1_example")
    genome <- PreprocessGenome("hg38", vcf = vcf_chr1_example)
    mut_table <- MutationsTable(vcf_chr1_example, genome, 3)
    type_table <- TypeCounts(mut_table)
    expect_equal(inherits(type_table, "table"), TRUE)
})
