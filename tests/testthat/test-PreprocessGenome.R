library(testthat)

test_that("Running PreprocessGenome with non-existent file", {
    path <- system.file("example-data", "not_exists.txt", package="SNVfromVCF")
    expect_error(PreprocessGenome(path), "The file does not exist.")
})

test_that("Running PreprocessGenome with a wrong file type",{
    path <- system.file("example-data", "test.txt", package="SNVfromVCF")
    expect_error(PreprocessGenome(path), "Reference genome file must be in fasta format")
})

test_that("Running PreprocessGenome with a .fasta file",{
    path <- system.file("example-data", "fakegenome.fasta", package="SNVfromVCF")
    genome <- PreprocessGenome(path)
    expect_equal(length(genome), 1)
    expect_equal(names(genome)[1], "chr1")
})

test_that("Running PreprocessGenome with a .fa file",{
    path <- system.file("example-data", "fakegenome.fa", package="SNVfromVCF")
    genome <- PreprocessGenome(path)
    expect_equal(length(genome), 1)
    expect_equal(names(genome)[1], "chr1")
})

test_that("Running PreprocessGenome with a .fna file",{
    path <- system.file("example-data", "fakegenome.fna", package="SNVfromVCF")
    genome <- PreprocessGenome(path)
    expect_equal(length(genome), 1)
    expect_equal(names(genome)[1], "chr1")
})

test_that("If chromosomes number is different from 1 or 25 after filtering", {
    path <- system.file("example-data", "fakegenome_with_2_chr.fna",
                        package="SNVfromVCF")
    expect_error(PreprocessGenome(path),
                 "Provide a full genome or a single chromosome")
})

test_that("Running PreprocessGenome with hg19", {
    genome <- PreprocessGenome("hg19")
    expect_equal(length(genome), 25)
    expect_equal(names(genome), c(paste0("chr", c(seq_len(22))),
                                  "chrX", "chrY", "chrM"))
})

test_that("Running PreprocessGenome with hg38", {
    genome <- PreprocessGenome("hg19")
    expect_equal(length(genome), 25)
    expect_equal(names(genome), c(paste0("chr", c(seq_len(22))),
                                  "chrX", "chrY", "chrM"))
})

test_that("Running PreprocessGenome filtering using a VCF object", {
    data("vcf_chr1_example")
    genome <- PreprocessGenome("hg38", vcf = vcf_chr1_example)
    expect_equal(length(genome), 1)
    expect_equal(names(genome), seqlevels(vcf_chr1_example))
})

