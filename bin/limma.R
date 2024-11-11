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
metavar = c('patient', 'cell_type', 'chemo', 'IO', 'TKI',  'V2V3',
            'age', 'sex', 'race', 'smoking', 'EGFR_status') # 'tissue',
obs_df = cts[,metavar]
obs_df$chemo = factor(1*(obs_df$chemo=='Chemo'), levels=c(0,1))
obs_df$IO = factor(1*(obs_df$IO=='IO'), levels=c(0,1))
obs_df$TKI = factor(1*(obs_df$TKI=='TKI'), levels=c(0,1))
#obs_df$LN = factor(1*(obs_df$tissue=='LN'), levels=c(0,1))
#obs_df$PDX = factor(1*(obs_df$tissue=='PDX'), levels=c(0,1))
#obs_df$metastasis = factor(1*(!obs_df$tissue %in% c('Lung','LN','PDX')), levels=c(0,1))
obs_df$V2V3 = factor(1*(obs_df$V2V3 == 'ten_x_v3'), levels=c(0,1))
obs_df$EGFR_status = factor(1*(obs_df$EGFR_status == 'Yes'), levels=c(0,1))
obs_df$sex = as.numeric(obs_df$sex == 'Male')
obs_df$age = as.numeric(obs_df$age)
obs_df$smoking = as.numeric(obs_df$smoking != 'never')
obs_df$race = as.numeric(obs_df$race !='White')
cell_types_split = strsplit(cell_types, "_")[[1]]
obs_df$cell_type = factor(match(obs_df$cell_type, cell_types_split), levels=1:length(cell_types_split))

cts = cts[, !(colnames(cts) %in% metavar)]
cts = t(cts[rowSums(cts) != 0, colSums(cts) != 0])

# Define regression formula and create design matrix
reg_formula = '~cell_type'
#V2V3 + LN + PDX + metastasis + chemo + IO + TKI + EGFR_status+ smoking + sex + race + 
design <- model.matrix(as.formula(reg_formula), obs_df)

# Filter and normalize data
dge <- DGEList(counts=cts)
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

# Perform voom transformation with block correlation
v <- voom(dge,design, plot=T)
corfit <- duplicateCorrelation(v,design,block=obs_df$patient)
print(corfit$consensus)

v <- voom(dge,design, plot=T, block=obs_df$patient, correlation = corfit$consensus)
corfit <- duplicateCorrelation(v,design,block=obs_df$patient)

# Fit linear model and get top table
fit <- lmFit(v, design, block=obs_df$patient, correlation=corfit$consensus)
fit <- eBayes(fit, trend=T, robust=T)
coef_names = colnames(coef(fit))[-1]
tt_filtered = data.frame()
for (x in coef_names) {
  tt = topTable(fit, coef=x, number=30000, sort.by='none')
  print(max(tt$logFC))
  print(min(tt$adj.P.Val))
  tt_filtered_temp = tt %>% 
    dplyr::filter(abs(logFC) > log2(1.5), adj.P.Val < 0.1) %>%
    dplyr::arrange(-logFC)
  tt_filtered_temp$coef = x
  tt_filtered = rbind(tt_filtered, tt_filtered_temp)
}

# Write results
write.table(tt, paste0(odir, '/limma.paired_', cell_types,'.txt'), sep='\t', quote=F, row.names=T, col.names=NA)
write.table(tt_filtered, paste0(odir, '/limma.paired_', cell_types, '.filtered.txt'), sep='\t', quote=F, row.names=T, col.names=NA)
