# Install ACE package and run ACE analysis
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("QDNEseq")
BiocManager::install("ACE")
BiocManager::install("QDNAseq.hg19")
install.packages("devtools")
devtools::install_github("asntech/QDNAseq.hg38@main")

# Setting
## specify the directory containing your bam-files
userpath <- "/Users/chun/WTL Dropbox/YenChun Chen/TumorEvolutionLab/Project_Peritoneal_Metastasis/wgs/O36/bam.GRCh38"
outputdir <- "/Users/chun/Documents/GitHub/ace-pipeline/results/OM1"
ref_genome <- "hg38"
bin <- 50

# Load libraries
library(ACE)
library(Biobase)
library(QDNAseq)
library(ggplot2)
library(BiocGenerics)
library(BiocParallel)
if (ref_genome == "hg19") {
  library(QDNAseq.hg19)
} else if (ref_genome == "hg38") {
  library(QDNAseq.hg38)
}

runACE(userpath,
  outputdir,
  filetype = "bam", binsizes = c(bin),
  ploidies = c(2, 3, 4, 5, 6), imagetype = "png", genome = ref_genome
)
