#
# Author: Mack Yuen
# Date: 2025-01-27
# Create PCA plot
#

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tibble)


# Read CPM file
logCPM <- read_delim("Results/Fraser_Fir_MeJa_logCPM_20250219.txt", 
                     col_names = TRUE, show_col_types = FALSE)

# Remove first column header to run PCA
logCPM <- as.data.frame(logCPM) |> 
  column_to_rownames('contig')

needles <- str_detect(colnames(logCPM), 'N')

needles_logCPM <- logCPM[needles]

pca <- prcomp(t(needles_logCPM), scale. = TRUE)

summary <- summary(pca)

pc1_var <- summary$importance[2, 1] * 100
pc2_var <- summary$importance[2, 2] * 100


dat <- as.data.frame(pca$x)

dat$genotype <- str_extract(rownames(dat), "G\\d\\d") |>
  as.factor()

dat$treatment <- str_split_i(rownames(dat), '', 4) |>
  as.factor()

dat$tissue <- str_split_i(rownames(dat), '', 5) |>
  as.factor()

dat$treatment <- str_split_i(rownames(dat), '', 4)

g <-  ggplot(dat, aes(PC1, PC2, label = rownames(dat))) +
  geom_text(vjust = 2.25, alpha = 0.25, size = 3) + 
  geom_point(aes(color = genotype, shape = treatment), size = 3) +
  labs(title = "Abies fraseri needles PCA",
       x = str_c("PC1 (", pc1_var, "%)"),
       y = str_c("PC2 (", pc2_var, "%)")) + 
  theme_bw()

plotname <- str_c("Results/Fraser_Fir_MeJa_needles_PCA.png")

ggsave(plotname, g)
