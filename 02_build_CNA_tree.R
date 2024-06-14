library(ACE)
library(ggplot2)
library(gplots)
library(ape)
library(phangorn)
library(RColorBrewer)

setwd("/Users/chun/Documents/GitHub/ace-pipeline/results/OM1")
dir.create("segmentdf")

object <- readRDS("50kbp.rds")
# define some variables to hold data
samples <- object@phenoData@data$name
# read info file
## set purity manually after reviewing ACE suggestions
## rationale: *clonal* alterations 1p loss and 4 loss are exactly at copy number 1
info <- read.delim("info.txt", sep = "\t")

samples_paper_name <- info$samples
samples <- samples_paper_name
patient <- "OM1"
segments <- list()
cellularities <- list()
bins <- list()
setpurity <- info$purity
# loop through data object and collect relevant bin-level information
ploidy <- 2
for (i in 1:length(samples)) {
  template <- objectsampletotemplate(object, index = i)
  model1 <- singlemodel(object, QDNAseqobjectsample = i, ploidy = ploidy)
  bestfit1 <- model1$minima[tail(which(model1$rerror == min(model1$rerror)), 1)]

  model <- singlemodel(template, ploidy = ploidy)
  singleplot(template, ploidy = ploidy, cellularity = setpurity[i], title = samples[i])
  ggsave(filename = paste(patient, "_auto_set_", samples[i], ".pdf", sep = ""))
  segmentdf <- getadjustedsegments(template, ploidy = ploidy, cellularity = setpurity[i])

  # get bin data
  bindata <- template
  bindata <- bindata[!is.na(bindata$copynumbers), ]
  bindata$chr <- as.character(bindata$chr)
  # filter out very small events
  bindata$intcopy <- NA
  bindata$meancopy <- NA

  # generate adjusted bin data

  for (j in 1:dim(segmentdf)[1]) {
    chr <- as.character(segmentdf$Chromosome[j])
    start <- segmentdf$Start[j]
    end <- segmentdf$End[j]

    bindata$intcopy[which(bindata$chr == chr & bindata$start >= start & bindata$end <= end)] <- segmentdf$Copies[j]
    bindata$meancopy[which(bindata$chr == chr & bindata$start >= start & bindata$end <= end)] <- segmentdf$Segment_Mean[j]
  }

  # store segmentdf as tsv
  write.table(segmentdf, file = paste0("segmentdf/", samples[i], "_segmentdf.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # generate bindata with -1 for deletion, 0 for normal ploidy, 1 for amplification
  bindata$simple <- NA
  bindata$simple[bindata$intcopy == ploidy] <- 0
  bindata$simple[bindata$intcopy < ploidy] <- -1
  bindata$simple[bindata$intcopy > ploidy] <- 1

  segments[[i]] <- segmentdf
  cellularities[i] <- model$minima[which.min(model$rerror)]
  bins[[i]] <- bindata
}

# collect data across all samples in a matrix format

s <- array(0, dim = c(dim(bins[[1]])[1], length(samples)))
for (i in 1:length(samples)) s[, i] <- bins[[i]]$intcopy
colnames(s) <- do.call(rbind, strsplit(samples, split = "_"))[, 1]

# make a heatmap

# c <- brewer.pal(8, "PiYG")

# pdf(paste0("copynumber-O36-full-manualset_heatmap2.pdf"), width = 10, height = 10)
# # draw(heat, ht_gap = unit(0, "cm"))
# heatmap.2(s, trace = "none", col = c, dendrogram = c("column"), Rowv = F, breaks = seq(0, 4, 0.5), key.title = NA, key.xlab = "copy number", lhei = c(2, 7), key.ylab = NA)
# dev.off()

# define useful tip colors make a NJ tree
colors <- samples

# MA2 colors
colors[grep("L1", colors)] <- "olivedrab4"
colors[grep("Col", colors)] <- "cyan4"
colors[grep("P", colors)] <- "firebrick3"
colors[grep("N", colors)] <- "black"
colors[grep("Ome", colors)] <- "orange2"
colors[grep("Stu", colors)] <- "sienna"
colors[grep("X", colors)] <- "skyblue4"
colors[grep("Bla", colors)] <- "royalblue4"
colors[grep("Bre", colors)] <- "palevioletred"

pdf("MA2_manualset_100kb_Euclidian_NJ.pdf")
d <- dist(t(s), method = "euclidian")
plot.phylo(nj(d), type = "unrooted", tip.color = colors)
dev.off()

# mean data for comparison
sm <- array(0, dim = c(dim(bins[[1]])[1], length(samples)))
for (i in 1:length(samples)) sm[, i] <- bins[[i]]$meancopy
colnames(sm) <- do.call(rbind, strsplit(samples, split = "_"))[, 1]

pdf("MA2_manualset_100kb_Euclidian_UPGMA_meancopy.pdf")
d <- dist(t(sm), method = "euclidian")
plot.phylo(upgma(d), type = "unrooted", tip.color = colors)
dev.off()


# Normal sample
normal_sample <- "N1"
source("utils/wGII_calculation_ace.R")
