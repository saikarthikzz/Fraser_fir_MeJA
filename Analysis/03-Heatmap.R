#
# Author: Mack Yuen
# Date: 2025-02-18
# Create Heatmap
#


library(heatmaply)
library(readr)
library(stringr)
library(dplyr)


# Read CPM file
logCPM <- read_delim("Results/Fraser_Fir_MeJa_logCPM_20250219.txt", 
                     col_names = TRUE, show_col_types = FALSE)

g51_de <- read_delim("Results/Fraser_Fir_g51_20250219.txt", 
                     col_names = TRUE, show_col_types = FALSE)

g85_de <- read_delim("Results/Fraser_Fir_g85_20250219.txt", 
                     col_names = TRUE, show_col_types = FALSE)

g97_de <- read_delim("Results/Fraser_Fir_g97_20250219.txt", 
                     col_names = TRUE, show_col_types = FALSE)

de_contigs <- c(g51_de$contig, g85_de$contig, g97_de$contig) |> 
  unique()

length(de_contigs)
# [1] 2246

needles <- str_detect(colnames(logCPM), 'N')

# We need to keep the contig names so we negate the first value
needles[1] <- TRUE

needles_logCPM <- logCPM[needles]

de_logCPM <- needles_logCPM |> 
  filter(contig %in% de_contigs)

dim(de_logCPM)
# [1] 2246   45

heatmaply(de_logCPM)
