## -----------------------------------------------------------------------------
library(EthSEQ)

## Run the analysis
ethseq.Analysis(
  target.vcf = system.file("extdata", "Samples.HGDP.10000SNPs.vcf",package="EthSEQ"),
  model.gds = system.file("extdata", "Reference.Gencode.Exome.10000SNPs.gds",package="EthSEQ"),
  out.dir = file.path(tempdir(),"EthSEQ_Analysis/"),
  verbose=TRUE,
  cores =1,
  composite.model.call.rate = 1,
  space = "3D")

## Load and display computed ethnicity annotations
ethseq.annotations = read.delim(file.path(tempdir(),"EthSEQ_Analysis/Report.txt"),
	sep="\t",as.is=TRUE,header=TRUE)
head(ethseq.annotations)

## Delete analysis folder
unlink(file.path(tempdir(),"EthSEQ_Analysis/"),recursive=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  library(EthSEQ)
#  
#  ## View all available reference models
#  getModelsList()
#  
#  ## Run the analysis
#  ethseq.Analysis(
#    target.vcf = system.file("extdata", "Samples.HGDP.10000SNPs.vcf",package="EthSEQ"),
#    model.available = "Gencode.Exome",
#    model.assembly = "hg38",
#    model.pop = "All",
#    out.dir = file.path(tempdir(),"EthSEQ_Analysis/"),
#    verbose=TRUE,
#    cores =1,
#    composite.model.call.rate = 1,
#    space = "3D")
#  
#  ## Delete analysis folder
#  unlink(file.path(tempdir(),"EthSEQ_Analysis/"),recursive=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  library(EthSEQ)
#  
#  ## Download BAM file used in the analysis
#  download.file("https://github.com/cibiobcg/EthSEQ_Data/raw/master/BAM/HGDP00228.sub_GRCh38.bam",
#                destfile = file.path(tempdir(),"HGDP00228.sub_GRCh38.bam"))
#  download.file("https://github.com/cibiobcg/EthSEQ_Data/raw/master/BAM/HGDP00228.sub_GRCh38.bam.bai",
#                destfile = file.path(tempdir(),"HGDP00228.sub_GRCh38.bam.bai"))
#  download.file("https://github.com/cibiobcg/EthSEQ_Data/raw/master/BAM/HGDP01200.sub_GRCh38.bam",
#                destfile = file.path(tempdir(),"HGDP01200.sub_GRCh38.bam"))
#  download.file("https://github.com/cibiobcg/EthSEQ_Data/raw/master/BAM/HGDP01200.sub_GRCh38.bam.bai",
#                destfile = file.path(tempdir(),"HGDP01200.sub_GRCh38.bam.bai"))
#  download.file("https://github.com/cibiobcg/EthSEQ_Data/raw/master/BAM/HGDP01201.sub_GRCh38.bam",
#                destfile = file.path(tempdir(),"HGDP01201.sub_GRCh38.bam"))
#  download.file("https://github.com/cibiobcg/EthSEQ_Data/raw/master/BAM/HGDP01201.sub_GRCh38.bam.bai",
#                destfile = file.path(tempdir(),"HGDP01201.sub_GRCh38.bam.bai"))
#  
#  ## Create BAM files list
#  write(c(file.path(tempdir(),"HGDP00228.sub_GRCh38.bam"),
#          file.path(tempdir(),"HGDP01200.sub_GRCh38.bam"),
#          file.path(tempdir(),"HGDP01201.sub_GRCh38.bam")),
#          file.path(tempdir(),"BAMs_List.txt"))
#  
#  ## Run the analysis
#  ethseq.Analysis(
#    bam.list = file.path(tempdir(),"BAMs_List.txt"),
#    model.available = "Gencode.Exome",
#    out.dir = file.path(tempdir(),"EthSEQ_Analysis/"),
#    verbose = TRUE,
#    cores = 1,
#    aseq.path = file.path(tempdir(),"EthSEQ_Analysis/"),
#    run.genotype = TRUE,
#    mbq = 20,
#    mrq = 20,
#    mdc = 10,
#    composite.model.call.rate = 1,
#    space = "3D",
#    bam.chr.encoding = TRUE) # chromosome names encoded without "chr" prefix in BAM files
#  
#  ## Delete analysis folder
#  unlink(file.path(tempdir(),"EthSEQ_Analysis/"),recursive=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  library(EthSEQ)
#  
#  ## Create multi-step refinement matrix
#  m = matrix("",ncol=2,nrow=2)
#  m[1,1] = "EUR|AFR|AMR"
#  m[2,2] = "EUR|AMR"
#  
#  ## Run the analysis on a toy example with only 10000 SNPs
#  ethseq.Analysis(
#    target.vcf = system.file("extdata","Samples.HGDP.10000SNPs.vcf",package="EthSEQ"),
#    out.dir = file.path(tempdir(),"EthSEQ_Analysis/"),
#    model.available = "Gencode.Exome",
#    verbose = TRUE,
#    refinement.analysis = m,
#    composite.model.call.rate = 1,
#    space = "3D")
#  
#  ## Delete analysis folder
#  unlink(file.path(tempdir(),"EthSEQ_Analysis/"),recursive=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  library(EthSEQ)
#  
#  ### Load list of VCF files paths
#  vcf.files = c(system.file("extdata","RefSample1.vcf", package="EthSEQ"),
#               	system.file("extdata","RefSample2.vcf", package="EthSEQ"))
#  
#  ### Load samples annotations
#  annot.samples = read.delim(system.file("extdata","Annotations_Test_v3.txt",package="EthSEQ"))
#  
#  ### Create reference model
#  ethseq.RM(
#    vcf.fn = vcf.files,
#    annotations = annot.samples,
#    out.dir = file.path(tempdir(),"EthSEQ_Analysis/"),
#    model.name = "Reference.Model",
#    bed.fn = NA,
#    call.rate = 1,
#    cores = 1)
#  
#  ## Delete example file
#  unlink(file.path(tempdir(),"EthSEQ_Analysis/"),recursive=TRUE)

