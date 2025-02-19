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
quant_files <- dir(path = 'Data',
                           pattern = 'quant.sf',
                           full.names = TRUE,
                           recursive = TRUE,
                           include.dirs = FALSE)

length(quant_files)
# [1] 44

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

tissues <- str_sub(samples, 5, 5) |> 
  as.factor()

biol_reps <- str_sub(samples, 6, 6)

exp_design <- data.frame(samples, genotypes, treatments, tissues, biol_reps)

exp_design$groups <- with(exp_design, interaction(genotypes, treatments))

exp_design
#    samples genotypes treatments tissues biol_reps groups
# 1   G51CB1       G51          C       B         1  G51.C
# 2   G51CB2       G51          C       B         2  G51.C
# 3   G51CB3       G51          C       B         3  G51.C
# 4   G51CN1       G51          C       N         1  G51.C
# 5   G51CN2       G51          C       N         2  G51.C
# 6   G51CN3       G51          C       N         3  G51.C
# 7   G51MB1       G51          M       B         1  G51.M
# 8   G51MB2       G51          M       B         2  G51.M
# 9   G51MB3       G51          M       B         3  G51.M
# 10  G51MN1       G51          M       N         1  G51.M
# You get the rest ...


# Create DGEList object
x <- DGEList(counts = expressions$counts, group = exp_design$groups)

colnames(x$counts) <- samples

# Remove contigs with low-expression
keepers <- filterByExpr(x, group = exp_design$groups)
  
table(keepers)
#  FALSE   TRUE
# 776320 156312
  
x <- x[keepers, ]

# Create model matrix
mod_matrix <- model.matrix(~ 0 + groups, exp_design)
colnames(mod_matrix) <- str_replace(colnames(mod_matrix), '\\.', '_')

mod_matrix
#    groupsG51_C groupsG85_C groupsG97_C groupsG51_M groupsG85_M groupsG97_M
# 1            1           0           0           0           0           0
# 2            1           0           0           0           0           0
# 3            1           0           0           0           0           0
# 4            1           0           0           0           0           0
# 5            1           0           0           0           0           0
# 6            1           0           0           0           0           0
# 7            0           0           0           1           0           0
# 8            0           0           0           1           0           0
# 9            0           0           0           1           0           0
# 10           0           0           0           1           0           0
# You get the rest ...


# TMM normalization
y <- calcNormFactors(x)
  

contrasts_matrix <- makeContrasts(
  g97 = groupsG97_M - groupsG97_C,
  g85 = groupsG85_M - groupsG85_C,
  g51 = groupsG51_C - groupsG51_M,
  levels = mod_matrix
)


# Transformation function
v <- cpm(y, log = TRUE)
  
log_cpm_filename <- str_c('Results/Fraser_Fir_MeJa_logCPM_', today, '.txt')
  
log_cpm <- data.frame(v) |> 
  rownames_to_column(var = "contig")
  
write_tsv(log_cpm, log_cpm_filename, num_threads = 4)

  
fit <- lmFit(v, mod_matrix)
  
fit2 <- contrasts.fit(fit, contrasts = contrasts_matrix)
  
fit3 <- eBayes(fit2, trend=TRUE)
  
summary(decideTests(fit3, method = 'separate', 
                    adjust.method = 'BH', 
                    p.value = 0.05, lfc = 2))

#           g97    g85    g51
# Down       96    459     12
# NotSig 156109 154220 156296
# Up        107   1633      4

# To get the contrasts
contrasts <- colnames(fit3$contrasts)

for (c in contrasts) {
  results <- topTable(fit3, coef = c, number = Inf, adjust.method = 'BH',
                      sort.by = 'logFC', p.value = 0.05, lfc = 2) |>
    rownames_to_column(var = 'contig')

  result_filename <- str_c('Results/Fraser_Fir_', c, '_', today, '.txt')

  write_tsv(results, result_filename, num_threads = 4)
}
