library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

data_dir <- "segmentdf/"

data <- list.files(data_dir, pattern = "\\.tsv$", full.names = TRUE)

index_dict <- list()
for (file in data) {
    sample <- str_split(file, "_")[[1]][1]

    # if (sample == normal_sample) {
    #   next
    # }

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
wGII_df <- data.frame(sample = names(index_dict), wGII = unlist(index_dict))

# Sort the dataframe by sample in ascending order
wGII_df <- wGII_df[order(wGII_df$sample), ]

# Save the dataframe to a CSV file
write.csv(wGII_df, file = paste0(patient, "_wGII.csv"), row.names = FALSE)

# Plot the wGII values
wGII_df <- wGII_df %>% arrange(desc(wGII))

ggplot() +
    geom_bar(data = wGII_df, aes(x = sample, y = wGII), fill = "#4793AF", alpha = 0.8) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "top"
    ) +
    ggsave(paste0(patient, "_wGII_plot.pdf"), dpi = 300)
