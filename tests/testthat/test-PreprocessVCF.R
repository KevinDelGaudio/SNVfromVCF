library(testthat)

test_that("Running PreprocessVCF with non-existent file", {
    path <- system.file("example-data", "not_exists.txt", package="SNVfromVCF")
    expect_error(PreprocessVCF(path), "The file does not exist.")
})

test_that("Running PreprocessVCF with a wrong file type",{
    path <- system.file("example-data", "test.txt", package="SNVfromVCF")
    expect_error(PreprocessVCF(path), "File must be in VCF format")
})

test_that("Running PreprocessVCF with a valid VCF file", {
    path <- system.file("example-data", "example.vcf", package="SNVfromVCF")
    vcf <- PreprocessVCF(path)
    refv <- as.character(ref(vcf))
    altv <- as.character(alt(vcf))
    expect_equal(typeof(vcf), "S4")
    expect_equal(class(vcf)[1], "ExpandedVCF")
    expect_equal(length(refv), length(altv))
    expect_equal(nchar(refv), rep(1, length(refv)))
    expect_equal(nchar(altv), rep(1, length(altv)))
})
