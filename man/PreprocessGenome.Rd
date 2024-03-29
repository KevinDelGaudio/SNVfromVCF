% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SNVfromVCF.R
\name{PreprocessGenome}
\alias{PreprocessGenome}
\title{Preprocess a genome for future use}
\usage{
PreprocessGenome(genome, vcf = NULL)
}
\arguments{
\item{genome}{(Mandatory) Genome name if using stock genomes provided by the
package or path to a custom fasta file. From a custom fasta files you can
provide only full genomes or a single chromosome; note that the headers of
the fasta file must be in the UCSC format or in the form ">chromosome 1,"
(the space and the comma are important, don't forget them!), this can be also
applied if you want to use a non canonical chromosome by modifying its
header.}

\item{vcf}{(Optional) An ExpandedVCF object useful for filtering the genome.}
}
\value{
filtered_genome A DNAStringSet object containing the sequences of
canonical chromosomes.
}
\description{
This function takes a genome and performs the required preprocessing steps
like keeping only standard chromosomes and set them up for use in other
functions of this package.
The genome can be loaded from a custom fasta file specified by the user or
directly from this function. If provided with a valid ExpandedVCF object this
function can use the names of chromosomes in the VCF file for further
filtering of the genome.
}
\examples{

genome <- PreprocessGenome(genome = "hg38") #load the assembly hg38

}
