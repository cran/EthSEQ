% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ethseq_RM.R
\name{ethseq.RM}
\alias{ethseq.RM}
\title{Create Reference Model for Ancestry Analysis}
\usage{
ethseq.RM(
  vcf.fn,
  annotations,
  out.dir = "./",
  model.name = "Reference.Model",
  bed.fn = NA,
  verbose = TRUE,
  call.rate = 1,
  cores = 1
)
}
\arguments{
\item{vcf.fn}{vector of paths to genotype files in VCF format}

\item{annotations}{data.frame with mapping of all samples names, ancestries and gender}

\item{out.dir}{Path to output folder}

\item{model.name}{Name of the output model}

\item{bed.fn}{path to BED file with regions of interest}

\item{verbose}{Print detailed information}

\item{call.rate}{SNPs call rate cutoff for inclusion in the final reference model}

\item{cores}{How many parallel cores to use in the reference model generation}
}
\value{
Logical value indicating the success of the analysis
}
\description{
This function creates a GDS reference model that can be used to performe EthSEQ ancestry analysis
}
