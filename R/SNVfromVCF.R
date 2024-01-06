#' Preprocess a VCF file
#'
#' This function loads a VCF file and filters it keeping only Single Nucleotide
#' Variants, it extends the CollapsedVCF object in order to have in the ALT
#' field only one variant per line then it return the ExpandedVCF object for
#' further analysis.
#'
#' @param vcf (Mandatory) Path to input file.
#' @return filtered_vcf An ExpandedVCF object.
#' @examples
#'
#' path <- system.file("example-data", "example.vcf", package="SNVfromVCF")
#' vcf <- PreprocessVCF(path)
#'
#' @export
PreprocessVCF <- function(vcf){
    #check if the file exists
    if (!file.exists(vcf)) {
        stop("The file does not exist.")
    }

    #check if the file has the correct extension
    if(!endsWith(vcf, ".vcf")){
        stop("File must be in VCF format")
    }

    #load the VCF file
    vcf <- VariantAnnotation::readVcf(vcf)

    #flatten the vcf file
    vcf <- VariantAnnotation::expand(vcf)

    #keep only SNVs
    allowed_variants <- c("A", "C", "G", "T")
    filtered_vcf <- vcf[as.character(ref(vcf)) %in% allowed_variants]
    filtered_vcf <- filtered_vcf[as.character(alt(filtered_vcf)) %in%
                                     allowed_variants]

    #the chromosome notation must be in the form chrN instead of only N
    seqlevels(rowRanges(filtered_vcf)) <- ifelse(
        !grepl("chr", seqlevels(rowRanges(filtered_vcf))),
        paste0("chr", seqlevels(rowRanges(filtered_vcf)) ),
        seqlevels(rowRanges(filtered_vcf)))


    return(filtered_vcf)
}

#' Preprocess a genome for future use
#'
#' This function takes a genome and performs the required preprocessing steps
#' like keeping only standard chromosomes and set them up for use in other
#' functions of this package.
#' The genome can be loaded from a custom fasta file specified by the user or
#' directly from this function. If provided with a valid ExpandedVCF object this
#' function can use the names of chromosomes in the VCF file for further
#' filtering of the genome.
#'
#' @param genome (Mandatory) Genome name if using stock genomes provided by the
#' package or path to a custom fasta file. From a custom fasta files you can
#' provide only full genomes or a single chromosome; note that the headers of
#' the fasta file must be in the UCSC format or in the form ">chromosome 1,"
#' (the space and the comma are important, don't forget them!), this can be also
#' applied if you want to use a non canonical chromosome by modifying its
#' header.
#' @param vcf (Optional) An ExpandedVCF object useful for filtering the genome.
#' @return filtered_genome A DNAStringSet object containing the sequences of
#' canonical chromosomes.
#' @examples
#'
#' genome <- PreprocessGenome(genome = "hg38") #load the assembly hg38
#'
#' @export
PreprocessGenome <- function(genome, vcf = NULL){
    #since it only makes sense to call variants on the non-alternative version
    #of a chromosome we filter all the other versions and scaffolds
    if(genome == "hg19"){
        genome <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
        filtered_genome <- BSgenome::getSeq(genome,
                        names=c(paste0("chr", c(seq_len(22))),
                                "chrX", "chrY", "chrM"))
    }
    else if(genome == "hg38"){
        genome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
        filtered_genome <- BSgenome::getSeq(genome,
                        names=c(paste0("chr", c(seq_len(22))),
                                "chrX", "chrY", "chrM"))
    }
    else{
        #check if the file exists
        if (!file.exists(genome)) {
            stop("The file does not exist.")
        }

        #check if it has the correct extension
        if(endsWith(genome, ".fasta") ||
           endsWith(genome, ".fna") ||
           endsWith(genome, ".fa")){
            genome <- Biostrings::readDNAStringSet(genome)
        }
        else{
            stop("Reference genome file must be in fasta format")
        }

        #filter to keep only non-alternative chromosomes
        filtered_genome <- genome[grepl(
            "(chr.+\\s([0-9]+|X|Y|MT)|mitochondrion),", names(genome))]

        #if given only one chromosome it preprocess it in order to have coherent
        #nomenclature of chromosomes, the same for full genomes from a custom
        #fasta file
        if(length(filtered_genome) == 1){
            if(grepl(".*chr.+\\s([0-9]+|X|Y|MT).*", names(filtered_genome)[1])){
            names(filtered_genome) <- gsub(".*chr.+\\s([0-9]+|X|Y|MT).*",
                                           "chr\\1", names(filtered_genome)[1])
            }
            else if(grepl(".*mitochondrion.*", names(filtered_genome)[1])){
                names(filtered_genome) <- gsub(".*mitochondrion.*",
                                               "chrM",names(filtered_genome)[1])
            }
        }
        else if(length(filtered_genome) == 25){
            names(filtered_genome) <- c(paste0("chr", c(seq_len(22))),
                                        "chrX", "chrY", "chrM")
        }
        else{
            stop("Provide a full genome or a single chromosome")
        }
    }

    #if a ExpandedVCF object is provided it uses the seqnames to filter further
    #the genome and not waste memory
    if (!is.null(vcf) && inherits(vcf, "ExpandedVCF")) {
        filtered_genome <- filtered_genome[seqnames(rowRanges(vcf))@values]
    }
    else if (!is.null(vcf) && !inherits(vcf, "ExpandedVCF")) {
        stop("Invalid ExpandedVCF object")
    }

    return(filtered_genome)
}

#' Compute a mutation table containing useful informations
#'
#' This function take as input a ExpandedVCF object, a genome and a context
#' length then it provides as output a dataframe containing for each mutation
#' chromosome, strand, position, type of SNV, the reference sequence and the
#' alternative sequence. Note that the mutation from G to A is the same as
#' C to T on the reverse strand, which implies a redundancy that is
#' addressed by converting the mutation types obtained for mutations of REF
#' bases A and G to their respective reverse complements.
#' Thus all mutation types are reported such that they have C or T as the
#' mutated REF base. The mutations which are kept as they are are on the
#' positive strand and the mutations for which the reverse complement is
#' computed are on the negative strand.
#'
#' @param vcf (Mandatory) ExpandedVCF object from which mutations are taken
#' @param genome (Mandatory) Genome used to retrieve the sequences
#' @param context_length (Mandatory, default: 3) The length of output sequences
#' @return mut_table data.frame containing the results
#' @examples
#'
#' data(vcf_chr1_example)
#' genome <- PreprocessGenome("hg38", vcf = vcf_chr1_example)
#' mut_table <- MutationsTable(vcf_chr1_example, genome, context_length = 3)
#'
#' @export
MutationsTable <- function(vcf, genome, context_length = 3){
    context_length <- as.integer(context_length)
    #check if vcf is an ExpandedVCF object
    if(!inherits(vcf, "ExpandedVCF")){
        stop("Invalid ExpandedVCF object")
    }
    #Check if genome is a DNAString
    if(!inherits(genome, "DNAStringSet")){
        stop("Provide a DNAStringSet genome")
    }
    #check if context length is a positive odd number
    if(context_length < 0 || (context_length %% 2) == 0) {
        stop("Parameter context_length must be a positive odd integer")
    }

    UP_DOWN <- (context_length - 1)/2
    CENTER <- UP_DOWN+1

    gr <- GenomicRanges::GRanges(seqnames = seqnames(vcf),
                  ranges = IRanges::IRanges(start = start(ranges(vcf))-UP_DOWN,
                                   end = start(ranges(vcf)) + UP_DOWN),
                  seqinfo = seqinfo(genome))

    ref <- getSeq(genome, gr)

    alt_df <- data.frame(old = as.character(ref), new = as.character(alt(vcf)))

    replace_sequence <- function(ref, alt, CENTER) {
        ref <- Biostrings::DNAString(ref)
        alt <- Biostrings::DNAString(alt)
        result <- as.character(Biostrings::replaceAt(ref,
                                IRanges::IRanges(CENTER, CENTER), alt))

        return(result)
    }

    alt_seq <- alt_df %>%
        rowwise() %>%
        mutate(alt_seq = replace_sequence(old, new, CENTER)) %>%
        ungroup()

    alt <- Biostrings::DNAStringSet(alt_seq$alt_seq)

    strand <- c()
    types <- c()
    for (i in seq_along(ref)){
        if(ref[[i]][CENTER] == Biostrings::DNAString("A") ||
           ref[[i]][CENTER] == Biostrings::DNAString("G")){
            ref[i] <- Biostrings::DNAStringSet(
                Biostrings::reverseComplement(ref[[i]]))
            alt[i] <- Biostrings::DNAStringSet(
                Biostrings::reverseComplement(alt[[i]]))
            strand <- append(strand, "-")
        }
        else{
            strand <- append(strand, "+")
        }
        refc <- as.character(ref[[i]][CENTER])
        altc <- as.character(alt[[i]][CENTER])
        types <- append(types, paste0("[", refc, ">", altc, "]"))
    }

    ref <- as.character(ref)
    type_df <- data.frame(first = substr(as.character(ref), 1, UP_DOWN),
                          second = substr(as.character(ref),
                                          CENTER+1,
                                          context_length),
                          type = types)

    write_type <- function(first, second, type){
        new_type <- paste0(first, type, second)
        return(new_type)
    }

    type_df2 <- type_df %>%
        rowwise() %>%
        mutate(new_type = write_type(first, second, type)) %>%
        ungroup()

    mut_table <- data.frame(
        CHR = as.character(seqnames(gr)),
        Strand = strand,
        Coordinates = as.character(ranges(gr)),
        MutationType = type_df2$new_type,
        REF = ref,
        ALT = as.character(alt)
    )

    return(mut_table)
}

#' Count the occurrence of each SNV and plot it.
#'
#' This function counts the occurrences of each SNV and, if requested by the
#' user, it shows a barplot summarizing the counts.
#'
#' @param mut_table (Mandatory) A dataframe generated by the MutationsTable
#' function
#' @param plot (Mandatory, default: FALSE) A bool that computes the plot if true
#' @return count_table A table object with the counts of the mutation types.
#' @examples
#'
#' data(vcf_chr1_example)
#' genome <- PreprocessGenome("hg38", vcf = vcf_chr1_example)
#' mut_table <- MutationsTable(vcf_chr1_example, genome, context_length = 3)
#' count_table <- TypeCounts(mut_table, plot = TRUE)
#' #compute count table and show the barplot
#'
#' @export
TypeCounts <- function(mut_table, plot = FALSE){
    if(is.data.frame(mut_table) &&
        ("CHR" %in% colnames(mut_table)) &&
        ("Strand" %in% colnames(mut_table)) &&
        ("Coordinates" %in% colnames(mut_table)) &&
        ("MutationType" %in% colnames(mut_table)) &&
        ("REF" %in% colnames(mut_table)) &&
        ("ALT" %in% colnames(mut_table))) {
            count_table <- table(mut_table$ALT)
    }
    else{
        stop("Invalid Mutations Table")
    }

    if(plot){
        df <- as.data.frame(count_table)
        p <- ggplot(df, aes(x = Var1, y = Freq, fill = Var1)) +
                        geom_bar(stat = "identity") +
                        scale_fill_hue(c = 40) +
                        labs(x = "Mutation Type", y = "Counts") +
                        ggtitle("Barplot of Mutation Types Counts") +
             theme(plot.title = element_text(hjust = 0.5),
                   legend.position = "none",
                   axis.text.x = element_text(angle = 90,
                                              vjust = 0.5,
                                              hjust=1))
        print(p)
    }

    return(count_table)
}
