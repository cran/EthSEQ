## ------------------------------------------------------------------------
library(EthSEQ)

out.dir = file.path(tempdir(),"EthSEQ_Analysis/")

## Run the analysis
ethseq.Analysis(
  target.vcf = system.file("extdata", "Samples_SS2_10000SNPs.vcf",package="EthSEQ"),
  out.dir = out.dir,
  model.gds = system.file("extdata", "Reference_SS2_10000SNPs.gds",package="EthSEQ"),
  verbose=TRUE,
  composite.model.call.rate = 1,
  space = "3D") # Default space is 2D

# Load and display computed ethnicity annotations
ethseq.annotations = read.delim(file.path(out.dir,"Report.txt"),sep="\t",as.is=TRUE,header=TRUE)
head(ethseq.annotations)
 
## Delete analysis folder
unlink(out.dir,recursive=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  library(EthSEQ)
#  
#  data.dir = file.path(tempdir(),"EthSEQ_Data/")
#  out.dir = file.path(tempdir(),"EthSEQ_Analysis/")
#  
#  ## Download genotype data in VCF format
#  dir.create(data.dir)
#  download.file("https://github.com/aromanel/EthSEQ_Data/raw/master/Sample_SS2.vcf",
#    destfile = file.path(data.dir,"Sample_SS2.vcf"))
#  
#  ## Run the analysis
#  ethseq.Analysis(
#    target.vcf =  file.path(data.dir,"Sample_SS2.vcf"),
#    out.dir = out.dir,
#    model.available = "SS2.Major",
#    model.folder = data.dir,
#    verbose=TRUE,
#    composite.model.call.rate = 1,
#    space = "3D") # Default space is 2D
#  
#  ## Delete analysis folder
#  unlink(data.dir,recursive=TRUE)
#  unlink(out.dir,recursive=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  library(EthSEQ)
#  
#  data.dir = file.path(tempdir(),"EthSEQ_Data")
#  out.dir = file.path(tempdir(),"EthSEQ_Analysis")
#  
#  ## Download BAM file used in the analysis
#  dir.create(data.dir)
#  download.file(
#   "https://github.com/aromanel/EthSEQ_Data/raw/master/NA07357_only10000SNPs.bam",
#   destfile = file.path(data.dir,"Sample.bam"))
#  download.file(
#   "https://github.com/aromanel/EthSEQ_Data/raw/master/NA07357_only10000SNPs.bam.bai",
#   destfile = file.path(data.dir,"Sample.bam.bai"))
#  
#  ## Create BAM files list
#  write(file.path(data.dir,"Sample.bam"),file.path(data.dir,"BAMs_List.txt"))
#  
#  ## Run the analysis
#  ethseq.Analysis(
#    bam.list = file.path(data.dir,"BAMs_List.txt"),
#    out.dir = out.dir,
#    model.gds = system.file("extdata","Reference_SS2_10000SNPs.gds",
#       package="EthSEQ"),
#    verbose=TRUE,
#    aseq.path = out.dir,
#    mbq=20,
#    mrq=20,
#    mdc=10,
#    run.genotype = TRUE,
#    composite.model.call.rate = 1,
#    cores=1,
#    bam.chr.encoding = FALSE) # chromosome names encoded without "chr" prefix in BAM files
#  
#  ## Delete analysis folder
#  unlink(data.dir,recursive=TRUE)
#  unlink(out.dir,recursive=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  library(EthSEQ)
#  
#  out.dir = file.path(tempdir(),"EthSEQ_Analysis")
#  data.dir = file.path(tempdir(),"EthSEQ_Data")
#  
#  ## Download genotype data in VCF format
#  dir.create(data.dir)
#  download.file("https://github.com/aromanel/EthSEQ_Data/raw/master/Target_SS2_10000SNPs.gds",
#    destfile = file.path(data.dir,"Target_SS2_10000SNPs.gds"))
#  
#  ## Create multi-step refinement matrix
#  m = matrix("",ncol=2,nrow=2)
#  m[1,1] = "SAS|EUR|EAS"
#  m[2,2] = "SAS|EUR"
#  
#  ## Run the analysis on a toy example with only 10000 SNPs
#  ethseq.Analysis(
#    target.gds = file.path(data.dir,"Target_SS2_10000SNPs.gds"),
#    out.dir = out.dir,
#    model.gds = system.file("extdata","Reference_SS2_10000SNPs.gds",
#  	package="EthSEQ"),
#    verbose=TRUE,
#    composite.model.call.rate = 1,
#    refinement.analysis = m,
#    space="3D")
#  
#  ## Delete analysis folder
#  unlink(out.dir,recursive=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  library(EthSEQ)
#  
#  out.dir = tempdir()
#  dir.create(out.dir)
#  
#  ### Load list of VCF files paths
#  vcf.files =
#    c(system.file("extdata", "VCF_Test_1.vcf", package="EthSEQ"),
#      system.file("extdata", "VCF_Test_2.vcf", package="EthSEQ"))
#  
#  ### Load samples annotations
#  annot.samples = read.delim(system.file("extdata", "Annotations_Test.txt",
#  	package="EthSEQ"))
#  
#  ### Create reference model
#  ethseq.RM(
#    vcf.fn = vcf.files,
#    annotations = annot.samples,
#    out.dir = out.dir,
#    model.name = "Reference.Model",
#    bed.fn = NA,
#    call.rate = 1,
#    cores = 1)
#  
#  ## Delete example file
#  unlink(out.dir,recursive=TRUE)

