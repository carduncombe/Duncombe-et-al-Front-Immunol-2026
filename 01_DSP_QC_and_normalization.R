# ---
# title: "GeoMx DSP Processing: QC, Normalization, and Batch Correction"
# author: "Caroline Duncombe"
# date: "2026-01-06"
# ---
  
  ## Overview
  
#   This script processes NanoString GeoMx DSP Whole Transcriptome Atlas (WTA) data
# from raw DCC files through quality control, normalization, ROI filtering,
# and batch correction. The final output is a Seurat object suitable for
# downstream differential expression analysis.
# 
# Pipeline summary:
# 1. Read GeoMx DCC files and metadata
# 2. Apply NanoString-recommended probe- and ROI-level QC
# 3. Normalize using housekeeping genes (GeomxTools) and log-normalization (Seurat)
# 4. Filter ROIs by size (200 µm) and liver zonation
# 5. Perform batch correction using ComBat-seq (sva)

# ---
  
## 1. Load libraries and define file paths
  

library(here)
library(GeomxTools)
library(Seurat)
library(dplyr)
library(tibble)
library(sva)


set.seed(123)

# Define input paths
dcc_dir        <- here("raw_data", "dcc_files")
pkc_file       <- here("raw_data", "pkc_file", "Mm_R_NGS_WTA_v1.0.pkc")
annotation_xls <- here("raw_data", "meta_data", "annotation_data_final.xlsx")

## 2- Read GeoMx DSP data

dcc_files <- dir(dcc_dir, pattern = "\\.dcc$", full.names = TRUE)

geomx_data <- suppressWarnings(
  readNanoStringGeoMxSet(
    dccFiles              = dcc_files,
    pkcFiles              = pkc_file,
    phenoDataFile         = annotation_xls,
    phenoDataSheet        = "Master",
    phenoDataDccColName   = "SAMPLE_ID",
    protocolDataColNames  = c("SegmentDisplayName", "RawReads", "unique_ID")
  )
)

### 3. Quality control filtering
## 3.1 Probe-level QC

geomx_qc <- setBioProbeQCFlags(
  geomx_data,
  qcCutoffs = list(minProbeRatio     = 0.1,
                   percentFailGrubbs = 20),
  removeLocalOutliers = FALSE
)

## 3.2 ROI-level QC

geomx_qc <- setSegmentQCFlags(
  geomx_qc,
  qcCutoffs = list(
    minSegmentReads    = 1000,
    percentAligned     = 80,
    percentStitched    = 80,
    percentTrimmed     = 80,
    percentSaturation  = 50,
    minNegativeCount   = 0,
    maxNTCCount        = 1000,
    minArea            = 1800,
    minNuclei          = 2
  )
)

# Identify ROIs passing all QC criteria
qc_flags <- protocolData(geomx_qc)[["QCFlags"]]
pass_idx <- which(rowSums(qc_flags) == 0)

geomx_qc_passed <- geomx_qc[, pass_idx]
dim(geomx_qc_passed) 
#Features  Samples 
#20175      273 


### 4. Count aggregation and normalization

# 4.1 Aggregate counts and housekeeping-gene normalization

geomx_counts <- aggregateCounts(geomx_qc_passed)

geomx_counts <- normalize(geomx_counts, norm_method="hk",fromElt="exprs", toElt="hk_norm")

# 4.2 Convert to Seurat object

seurat_obj <- as.Seurat(geomx_counts, normData = "exprs", forceRaw = TRUE)

### 5. Log-normalization and dimensional reduction (Seurat)

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 5000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

### 6. ROI filtering by size and zonation

# Retain only 200 µm ROIs
seurat_200um <- subset(seurat_obj, ROI_size == "200um")

# Exclude ROIs without zonation assignment
seurat_200um <- subset(seurat_200um, Liver_zone != "exclude")

dim(seurat_200um)
table(seurat_200um@meta.data$Group, seurat_200um@meta.data$Liver_zone, seurat_200um@meta.data$Type)

### 7. Batch correction using ComBat-seq (sva)

# Load dataset
dat <- seurat_200um

# Extract raw count matrix
count_matrix <- GetAssayData(dat, assay = "GeoMx", layer = "counts")
count_matrix <- as.matrix(count_matrix)

stopifnot(all(colnames(count_matrix) == rownames(dat@meta.data)))

# Define batch and covariates
batch <- factor(dat@meta.data$Block)

covariates <- data.frame(
  sex       = dat@meta.data$Group,
  infection = dat@meta.data$Type
)

# Apply ComBat-seq
combat_counts <- ComBat_seq(
  count_matrix,
  batch     = batch,
  group     = NULL,
  covar_mod = covariates
)

# Store batch-corrected counts as new assay
dat[["SVA_batchCorrected"]] <- CreateAssayObject(
  counts = as.matrix(combat_counts)
)

svacounts <- GetAssayData(dat, assay = "SVA_batchCorrected", layer = "counts") #pulling from object
dim(svacounts) # Checking that matrix is there. 



### SAVE
saveRDS(dat,file = here("clean_data","srt","QC_normalized_data_200um_SVA.rds"))


###### ALL DONE!!! ######

