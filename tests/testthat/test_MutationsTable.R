library(testthat)

test_that("Invalid ExpandedVCF, valid Genome, odd context_length",{
    vcf <- c(1) #not and ExpandedVCF object
    genome <- PreprocessGenome("hg38")
    expect_error(MutationsTable(vcf, genome, 3),
                 "Invalid ExpandedVCF object")
})

test_that("Valid ExpandedVCF, invalid Genome, odd context_length",{
    data("vcf_chr1_example")
    genome <- c(1) #not a DNAStringSet
    expect_error(MutationsTable(vcf_chr1_example, genome, 3),
                 "Provide a DNAStringSet genome")
})

test_that("Valid ExpandedVCF, valid Genome, even context_length",{
    data("vcf_chr1_example")
    genome <- PreprocessGenome("hg38", vcf = vcf_chr1_example)
    context_length <- 2 #even number
    expect_error(MutationsTable(vcf_chr1_example, genome, context_length),
                 "Parameter context_length must be a positive odd integer")
})

test_that("Checking if the output dataframe is correct",{
    data("vcf_chr1_example")
    genome <- PreprocessGenome("hg38", vcf = vcf_chr1_example)
    mut_table <- MutationsTable(vcf_chr1_example, genome, context_length = 3)
    valid <- (is.data.frame(mut_table) &&
                ("CHR" %in% colnames(mut_table)) &&
                ("Strand" %in% colnames(mut_table)) &&
                ("Coordinates" %in% colnames(mut_table)) &&
                ("MutationType" %in% colnames(mut_table)) &&
                ("REF" %in% colnames(mut_table)) &&
                ("ALT" %in% colnames(mut_table)))
    expect_equal(valid, TRUE)
    expect_match(mut_table$MutationType, ".*\\[[CT]>[ACGT]\\].*", all = TRUE)
})
