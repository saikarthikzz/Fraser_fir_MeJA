#
# Author: Mack Yuen
# Date: 2025-01-27
# Differential expression analysis with limma
#

library(edgeR)
library(limma)
library(readr)
library(stringr)
library(tibble)
library(tximport)


# Nark all output by date
today <- Sys.Date()
today <- format(today, "%Y%m%d")


# Salmon quant file location
quant_files <- dir(path = 'Needle-Data',
                           pattern = 'quant.sf',
                           full.names = TRUE,
                           recursive = TRUE,
                           include.dirs = FALSE)

# Remove outliers from post DE analysis
quant_files <- quant_files[!str_detect(quant_files, '97CN3')]
quant_files <- quant_files[!str_detect(quant_files, '97MN3')]

quant_files <- quant_files[!str_detect(quant_files, '85CN3')]
quant_files <- quant_files[!str_detect(quant_files, '85MN3')]

length(quant_files)
# [1] 18


# Read expressions with tximport
expressions <- tximport(
  quant_files,
  type = "salmon",
  txIn = TRUE,
  txOut = TRUE,
  countsFromAbundance = "lengthScaledTPM")
  

# Setup experimental design. Extract information from file names to avoid typos.
samples <- str_split_i(quant_files, '/', 2)

samples <- str_split_i(samples, '_', 1)

# Using numbers as the first letter of a variable is not a good idea.
# Patch sample name with a letter. Will save a lot of hassle down the road.
samples <- str_c('G', samples)

genotypes <- str_sub(samples, 1, 3) |> 
  as.factor()

treatments <- str_sub(samples, 4, 4) |> 
  as.factor()

biol_reps <- str_sub(samples, 6, 6)

exp_design <- data.frame(samples, genotypes, treatments, biol_reps)

exp_design$groups <- with(exp_design, interaction(genotypes, treatments))

exp_design
#    samples genotypes treatments biol_reps groups
# 1   G51CN1       G51          C         1  G51.C
# 2   G51CN2       G51          C         2  G51.C
# 3   G51CN3       G51          C         3  G51.C
# 4   G51MN1       G51          M         1  G51.M
# 5   G51MN2       G51          M         2  G51.M
# 6   G51MN3       G51          M         3  G51.M
# 7   G85CN1       G85          C         1  G85.C
# 8   G85CN2       G85          C         2  G85.C
# 9   G85CN4       G85          C         4  G85.C
# 10  G85MN1       G85          M         1  G85.M
# You get the rest ...


# Create DGEList object
x <- DGEList(counts = expressions$counts, group = exp_design$groups)

colnames(x$counts) <- samples

# Remove contigs with low-expression
keepers <- filterByExpr(x, group = exp_design$groups)
  
table(keepers)
# FALSE  TRUE 
# 89387 57302  

x <- x[keepers, ]

# Create model matrix
mod_matrix <- model.matrix(~ 0 + groups, exp_design)
colnames(mod_matrix) <- str_replace(colnames(mod_matrix), '\\.', '_')

mod_matrix
#    groupsG51_C groupsG85_C groupsG97_C groupsG51_M groupsG85_M groupsG97_M
# 1            1           0           0           0           0           0
# 2            1           0           0           0           0           0
# 3            1           0           0           0           0           0
# 4            0           0           0           1           0           0
# 5            0           0           0           1           0           0
# 6            0           0           0           1           0           0
# 7            0           1           0           0           0           0
# 8            0           1           0           0           0           0
# 9            0           1           0           0           0           0
# 10           0           0           0           0           1           0
# You get the rest ...


# TMM normalization
y <- calcNormFactors(x)
  

contrasts_matrix <- makeContrasts(
  g97 = groupsG97_M - groupsG97_C,
  g85 = groupsG85_M - groupsG85_C,
  g51 = groupsG51_M - groupsG51_C,
  levels = mod_matrix
)


# Transformation function
v <- cpm(y, log = TRUE)
  
log_cpm_filename <- str_c('Results/Fraser_Fir_MeJa_Needles_logCPM_', today, '.txt')
  
log_cpm <- data.frame(v) |> 
  rownames_to_column(var = "contig")
  
write_tsv(log_cpm, log_cpm_filename, num_threads = 4)

  
fit <- lmFit(v, mod_matrix)
  
fit2 <- contrasts.fit(fit, contrasts = contrasts_matrix)
  
fit3 <- eBayes(fit2, trend=TRUE)
  
summary(decideTests(fit3, method = 'separate', 
                    adjust.method = 'BH', 
                    p.value = 0.05, lfc = 2))

#          g97   g85   g51
# Down     116   108    17
# NotSig 55980 56304 57172
# Up      1206   890   113


# To get the contrasts
contrasts <- colnames(fit3$contrasts)

for (c in contrasts) {
  results <- topTable(fit3, coef = c, number = Inf, adjust.method = 'BH',
                      sort.by = 'logFC', p.value = 0.05, lfc = 2) |>
    rownames_to_column(var = 'contig')

  result_filename <- str_c('Results/Fraser_Fir_', c, '_', today, '.txt')

  write_tsv(results, result_filename, num_threads = 4)
}
