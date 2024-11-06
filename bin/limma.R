# Load necessary libraries
library(dplyr)
library('Matrix')
library(tibble)
library(zoo)
library(stringr)
library(lubridate)
library(dplyr)
library(tibble)
library(rhdf5)
library(limma)
library(edgeR)
library(commandr)

# Define command line arguments
args <- commandArgs(trailingOnly = TRUE)
opt <- commandr::commandArgs(quiet = TRUE)

# Define input and output directories
input_dir <- opt$input_dir
input_name <- opt$input_name
cell_types <- opt$cell_types
odir <- opt$output_dir

# Read data
main_dir = input_dir
fn = paste0(main_dir, '/obs.', input_name, '.csv')
obs_df = read.table(fn, sep = ',', header = T, row.names=1)
fn = paste0(main_dir, 'var_names', input_name, '.csv')
g = readLines(fn)

fn = paste0(main_dir, 'obs_names.', intput_name, '.csv', )
bc = readLines(fn)

cts = readMM(paste0(main_dir, 'counts.' input_name, '.mtx'))

cts = as.matrix(cts)

colnames(cts) = g
rownames(cts) = bc

cts = t(cts)
cts = cts[rowSums(cts) != 0,]

obs_df = obs_df[colnames(cts),]

obs_df = obs_df %>% dplyr::filter(colSums(cts) != 0)
cts = cts[,colSums(cts) != 0]

obs_df$patient = factor(gsub('_.*','',rownames(obs_df)))
dir.create(odir)

obs_df$histo = factor(obs_df$histo, levels=cell_types)

# Define regression formula
reg_formula = '~histo'

# Create design matrix
design <- model.matrix(as.formula(reg_formula), obs_df)
dge <- DGEList(counts=cts)
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

# Perform voom transformation
v <- voom(dge,design, plot=T)
corfit <- duplicateCorrelation(v,design,block=obs_df$patient)
print(corfit$consensus)

# Perform voom transformation with block correlation
v <- voom(dge,design, plot=T, block=obs_df$patient, correlation = corfit$consensus)
corfit <- duplicateCorrelation(v,design,block=obs_df$patient)
print(corfit$consensus)

# Fit linear model
fit <- lmFit(v, design, block=obs_df$patient, correlation=corfit$consensus)
fit <- eBayes(fit, trend=T, robust=T)

# Get top table
tt = topTable(fit, number=30000, sort.by='none')
cell_types = paste(cell_types, collapse='_')
write.table(tt, paste0(odir, '/limma.paired_', cell_types,'.txt',), sep='\t', quote=F, row.names=T, col.names=NA)
write.table(tt %>% dplyr::filter(adj.P.Val < 0.1 & abs(logFC) > log2(1.5)), 
            paste0(odir, '/limma.paired_', cell_types, '.filtered.txt'), sep='\t', quote=F, row.names=T, col.names=NA)
