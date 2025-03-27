#
# Author: Mack Yuen
# Date: 2025-02-18
# Create Heatmap
#


library(heatmaply)
library(readr)
library(stringr)
library(dplyr)
library(tibble)
library(purrr)


# Read CPM file
logCPM <- read_delim("Results/Fraser_Fir_MeJa_Needles_logCPM_20250320.txt", 
                     col_names = TRUE, show_col_types = FALSE)

g51_de <- read_delim("Results/Fraser_Fir_g51_20250320.txt", 
                     col_names = TRUE, show_col_types = FALSE)

g85_de <- read_delim("Results/Fraser_Fir_g85_20250320.txt", 
                     col_names = TRUE, show_col_types = FALSE)

g97_de <- read_delim("Results/Fraser_Fir_g97_20250320.txt", 
                     col_names = TRUE, show_col_types = FALSE)


# Concatenate all DE contigs
de_contigs <- c(g51_de$contig, g85_de$contig, g97_de$contig) |> 
  unique()

# Check the number of unique contigs that were differentially expressed
# across all genotypes
length(de_contigs)
# [1] 2794

# Filter logCPM for those only differentially expressed
de_logCPM <- logCPM |> 
  filter(contig %in% de_contigs)

dim(de_logCPM)
# [1] 2794   23


# We will take the column name and extract the genotype and treatment
# Since this is a tibble, we don't want the contig column
treatment <- colnames(de_logCPM)[-1]
treatment <- str_sub(treatment, 1, 5) |> 
  unique()

# Calculate the mean of logCPM by treatment
mean_de_logCPM <- map(treatment, function(t){
  df <- de_logCPM |> 
    select(starts_with(t))
  rowMeans(df)
})

names(mean_de_logCPM) <- treatment

mean_de_logCPM <- as.data.frame(mean_de_logCPM)

# Add back rownames to data frame
rownames(mean_de_logCPM) <- de_logCPM$contig

dim(mean_de_logCPM)
# [1] 2794    6

today <- Sys.Date()
today <- format(today, "%Y%m%d")

heatmaply(mean_de_logCPM, 
          scale = "row", 
          show_dendrogram = c(TRUE, FALSE), 
          cexRow = FALSE, 
          showticklabels = c(TRUE, FALSE),
          file = str_c("Results/draft_heatmap_", today, ".png"))
