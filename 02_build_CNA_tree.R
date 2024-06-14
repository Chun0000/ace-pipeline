library(ACE)
library(ggplot2)
library(gplots)
library(ape)
library(phangorn)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

setwd("/Users/chun/WTL Dropbox/YenChun Chen/TumorEvolutionLab/Project_Peritoneal_Metastasis/wgs/O36/copy.number.tree.result/GRCh38/50kbp")
dir.create("segmentdf")
dir.create("diagrams")

object <- readRDS("50kbp.rds")
# define some variables to hold data
samples <- object@phenoData@data$name

# read info file
## set purity manually after reviewing ACE suggestions
## rationale: *clonal* alterations 1p loss and 4 loss are exactly at copy number 1
info <- read.delim("info.txt", sep = "\t")
ploidy <- 2
patient <- "O36"
normal_sample <- "N1"
cal_normal_wGII <- TRUE

samples_paper_name <- info$samples
samples <- samples_paper_name
segments <- list()
cellularities <- list()
bins <- list()
setpurity <- info$purity

# loop through data object and collect relevant bin-level information
for (i in 1:length(samples)) {
  template <- objectsampletotemplate(object, index = i)
  model1 <- singlemodel(object, QDNAseqobjectsample = i, ploidy = ploidy)
  bestfit1 <- model1$minima[tail(which(model1$rerror == min(model1$rerror)), 1)]

  model <- singlemodel(template, ploidy = ploidy)
  singleplot(template, ploidy = ploidy, cellularity = setpurity[i], title = samples[i])
  ggsave(filename = paste("diagrams/", patient, "_auto_set_", samples[i], ".pdf", sep = ""))
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
color_tip <- info$color_tip

pdf("MA2_manualset_100kb_Euclidian_NJ.pdf")
d <- dist(t(s), method = "euclidian")
plot.phylo(nj(d), type = "unrooted", tip.color = color_tip)
dev.off()

# mean data for comparison
sm <- array(0, dim = c(dim(bins[[1]])[1], length(samples)))
for (i in 1:length(samples)) sm[, i] <- bins[[i]]$meancopy
colnames(sm) <- do.call(rbind, strsplit(samples, split = "_"))[, 1]

pdf("MA2_manualset_100kb_Euclidian_UPGMA_meancopy.pdf")
d <- dist(t(sm), method = "euclidian")
plot.phylo(upgma(d), type = "unrooted", tip.color = color_tip)
dev.off()

# calculate wGII

data_dir <- "segmentdf"

data <- list.files(data_dir, pattern = "\\.tsv$", full.names = TRUE)

index_dict <- list()
for (file in data) {
  sample <- str_split(file, "_")[[1]][1]

  print(sample)
  if (!cal_normal_wGII) {
    if (sample == normal_sample) {
      next
    }
  }

  df <- read.delim(file, sep = "\t")

  # Calculate length of each segment
  df$length <- df$End - df$Start

  # Identify aberrant segments
  df$is_aberrant <- df$Copies != ploidy

  # Calculate the total length of aberrant segments per Chromosome
  aberrant_lengths <- df[df$is_aberrant, ] %>%
    group_by(Chromosome) %>%
    summarise(aberrant_length = sum(length))

  # Calculate the total length of each Chromosome
  total_lengths <- df %>%
    group_by(Chromosome) %>%
    summarise(total_length = sum(length))

  # Merge to get a complete dataframe
  merged_lengths <- total_lengths %>%
    left_join(aberrant_lengths, by = "Chromosome") %>%
    replace_na(list(aberrant_length = 0))

  # Calculate the percentage of aberrant length per Chromosome
  merged_lengths$percentage_aberrant <- merged_lengths$aberrant_length / merged_lengths$total_length

  # Calculate the mean percentage aberration (weighted GII)
  wGII <- mean(merged_lengths$percentage_aberrant)

  # Store the wGII value
  index_dict[[sample]] <- wGII
}

# Convert the list to a dataframe
wGII_df <- data.frame(samples = names(index_dict), wGII = unlist(index_dict))

# Sort the dataframe by sample in ascending order
wGII_df <- wGII_df[order(wGII_df$samples), ]

# Save the dataframe to a CSV file
write.csv(wGII_df, file = paste0(patient, "_wGII.csv"), row.names = FALSE)

# Plot the wGII values
merged_data <- merge(wGII_df, info, by = "samples")
merged_data <- merged_data %>% arrange(desc(wGII))

ggplot(data = merged_data, aes(x = samples, y = wGII, fill = color_bar)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    legend.position = "none" # Remove the legend
  ) +
  scale_x_discrete(limits = merged_data$samples, name = "Sample") +
  scale_fill_manual(values = merged_data$color_bar) +
  scale_y_continuous(expand = expansion(add = c(0, 0.05)))

ggsave(paste0(patient, "_wGII_plot.pdf"), dpi = 300)
