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


# We need to keep the contig names so we negate the first value
# needles[1] <- TRUE

# needles_logCPM <- logCPM[needles]

de_logCPM <- logCPM |> 
  filter(contig %in% de_contigs)

dim(de_logCPM)
# [1] 2794   23

heatmaply(de_logCPM, draw_cellnote = FALSE, Colv = FALSE, RowSideColors = FALSE)


# We will take the column name and extract the genotype and treatment
# Since this is a tibble, we don't want the contig column
# mDat <- head(de_logCPM, 500)
treatment <- colnames(de_logCPM)[-1]
treatment <- str_sub(treatment, 1, 5) |> 
  unique()


# de_logCPM <- de_logCPM |> 
#   mutate(G51CN = mean(c(G51CN1, G51CN2, G51CN3)), .keep = "unused")
# 
# de_logCPM <- de_logCPM |> 
#   mutate(G51MN = mean(c(G51MN1, G51MN2, G51MN3)), .keep = "unused")
# 
# de_logCPM <- de_logCPM |> 
#   mutate(G85CN = mean(c(G85CN1, G85CN2, G85CN3)), .keep = "unused")
# 
# de_logCPM <- de_logCPM |> 
#   mutate(G51CN = mean(c(G51CN1, G51CN2, G51CN3)), .keep = "unused")


tmp <- map(treatment, function(t){
  df <- de_logCPM |> 
    select(starts_with(t))
  rowMeans(df)
})

names(tmp) <- treatment

tmp <- as.data.frame(tmp)

dim(tmp)

heatmaply(tmp, scale = "row", show_dendrogram = c(TRUE, FALSE), cexRow = FALSE, 
          showticklabels = c(TRUE, FALSE),
          file = "draft_heatmap_20250324.png")
