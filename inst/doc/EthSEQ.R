## ------------------------------------------------------------------------
library(EthSEQ)

## Run the analysis
ethseq.Analysis(
  target.vcf = system.file("extdata", "Samples_SS2_10000SNPs.vcf",
	package="EthSEQ"),
  out.dir = "/tmp/EthSEQ_Analysis/",
  model.gds = system.file("extdata", "Reference_SS2_10000SNPs.gds",
	package="EthSEQ"),
  verbose=TRUE,
  composite.model.call.rate = 1)

## Load and display computed ethnicity annotations
ethseq.annotations = read.delim("/tmp/EthSEQ_Analysis/Report.txt",
	sep="\t",as.is=TRUE,header=TRUE)
head(ethseq.annotations)

## Delete analysis folder
unlink("/tmp/EthSEQ_Analysis/",recursive=TRUE)

