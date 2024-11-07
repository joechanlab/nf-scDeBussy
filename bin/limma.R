# Load necessary libraries
library(dplyr)
library(Matrix)
library(tibble)
library(zoo)
library(stringr)
library(lubridate)
library(rhdf5)
library(limma)
library(edgeR)

# Define command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Define input and output directories
input_dir <- args[1]
cell_types <- args[2]
odir <- args[3]
if (!dir.exists(odir)) {
  dir.create(odir)
}

# Read data
cts = read.table(paste0(input_dir, '/counts.', cell_types, '.csv'), sep = ',', header = T, row.names=1)
obs_df = cts[,c('sample','cell_type')]
cts = cts[, !(colnames(cts) %in% c('sample', 'cell_type'))]
cts = t(cts[rowSums(cts) != 0, colSums(cts) != 0])
obs_df$cell_type = factor(obs_df$cell_type, levels=strsplit(cell_types, "_")[[1]])

# Define regression formula and create design matrix
reg_formula = '~cell_type'
design <- model.matrix(as.formula(reg_formula), obs_df)

# Filter and normalize data
dge <- DGEList(counts=cts)
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

# Perform voom transformation with block correlation
v <- voom(dge,design, plot=T)
corfit <- duplicateCorrelation(v,design,block=obs_df$sample)
print(corfit$consensus)

v <- voom(dge,design, plot=T, block=obs_df$sample, correlation = corfit$consensus)
corfit <- duplicateCorrelation(v,design,block=obs_df$sample)

# Fit linear model and get top table
fit <- lmFit(v, design, block=obs_df$sample, correlation=corfit$consensus)
fit <- eBayes(fit, trend=T, robust=T)
tt = topTable(fit, number=30000, sort.by='none')

# Write results
write.table(tt, paste0(odir, '/limma.paired_', cell_types,'.txt'), sep='\t', quote=F, row.names=T, col.names=NA)
write.table(tt %>% dplyr::filter(adj.P.Val < 0.1 & abs(logFC) > log2(1.5)), 
            paste0(odir, '/limma.paired_', cell_types, '.filtered.txt'), sep='\t', quote=F, row.names=T, col.names=NA)
