library('TCGAbiolinks')
library('easybio')
library('data.table')
library('SummarizedExperiment')
library('DESeq2')
library('tidyverse')

# Following this vignette: 
# https://cran.r-project.org/web/packages/easybio/vignettes/example_limma.html

# Get query (metadata) for downloading
query <- GDCquery(
  project = "TCGA-CHOL",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification"
)
# ============= exploring metadata =============
metadata <- query[[1]][[1]]
# Just checking what sample types there are
unique(metadata$sample_type)
# There's only Primary Tumor and Solid Tissue Normal

# Experimental stratery (RNA-seq, microarray?...)
unique(metadata$experimental_strategy)
# All RNA-Seq. Wonder if people still use microarray here...

# Analysis workflow?
unique(metadata$analysis_workflow_type)
# They're all STAR - Counts. Makes sense that they use STAR aligner for RNA-seq
# Seemed like they don't use HISAT2 aligner or Kallisto/Salmon pseudoaligner...

# Get number of samples
length(unique(metadata$cases))
# 44 samples total

# Get number of patients
length(unique(metadata$cases.submitter_id))
# 36 patients


# ============= Download ==================
# Download count data?
GDCdownload(query = query)

# Prepare into an R object?
# Convert into a RangedSumarizeExperiment object
# How is this different from SummarizeExperiment?
data <- GDCprepare(query = query)

# Some exploration of count data
assay(data)

# extracts and processes the necessary information from the TCGA data object, 
# separating tumor and non-tumor samples.
lt <- prepare_tcga(data)

lt$all$sampleInfo[["group"]] <- fifelse(lt$all$sampleInfo$sample_type == "Primary Tumor", "Tumor", "Normal")

colData(vsd)$group <- lt$all$sampleInfo$group

# --------------
dds <- DESeqDataSetFromMatrix(
  countData = assay(data, "unstranded"),
  colData   = colData(data),
  design    = ~ 1
)

dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)

pca <- prcomp(t(assay(vsd)))

pca_df <- data.frame(
  sample = colnames(vsd),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  group = colData(vsd)$group
)

ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  labs(
    title = "PCA (VST normalized)",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2,2], 1), "%)")
  )
# ---------------




