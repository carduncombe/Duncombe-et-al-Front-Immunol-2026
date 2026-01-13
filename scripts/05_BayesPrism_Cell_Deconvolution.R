

#BayesPrism Cellular deconvolution

#Conducted on November 13th, 2024.
#CONDUCTED ON 200um only samples. 

library(Seurat)
library(BayesPrism)
library(InstaPrism)
#devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")
#devtools::install_github("humengying0907/InstaPrism")

rm(list = ls())

# LOAD IN DATASETS -----------------------------------------------------------------------------------------------------------------------------

#SC dataset from published manuscript on Murine Liver
RDS <- readRDS("scRNAseq_from_paper/final_merged_curated_annotations_270623.RDS") # This is the snRNAseq data from Hildebrant et al. manuscript

# Normalized data with no batch correction:
dat_SUV <- readRDS("clean_data/srt/QC_normalized_data_200um_SVA.rds")


# CREATE THE OPTIONS FOR DOWN STREAM BAYES -----------------------------------------------------------------------------------------------------------------------------

#SUV combined method for infection types
dat_SUV_200_schizont <- subset(dat_SUV, subset = Type == "schizont")
dat_SUV_200_mock <- subset(dat_SUV, subset = Type == "mock")
dat_SUV_200_bystander <- subset(dat_SUV, subset = Type == "bystander")

#Further look at contribution by zone. 
dat_SUV_200_schizont_PP <- subset(dat_SUV_200_schizont, subset = Liver_zone == "PP")
dat_SUV_200_schizont_IZ <- subset(dat_SUV_200_schizont, subset = Liver_zone == "IZ")
dat_SUV_200_schizont_PC <- subset(dat_SUV_200_schizont, subset = Liver_zone == "PC")

dat_to_run_SUV <- list(
  dat_SUV_200_schizont = dat_SUV_200_schizont,
  dat_SUV_200_mock = dat_SUV_200_mock,
  dat_SUV_200_bystander = dat_SUV_200_bystander,
  dat_SUV_200_schizont_PP = dat_SUV_200_schizont_PP,
  dat_SUV_200_schizont_PC = dat_SUV_200_schizont_PC,
  dat_SUV_200_schizont_IZ = dat_SUV_200_schizont_IZ,
  dat_SUV_200 = dat_SUV
  )

# Calculate the average per cell type from the single cell dataset from the paper  -----------------------------------------------------------------------------------------------------------------------------

# Single-cell data
sc.dat <- AverageExpression(RDS, group.by = "merged_type") 
sc.dat <- t(as.matrix(sc.dat$RNA))

# Labels
cell.type.labels <- rownames(sc.dat)
cell.state.labels <- rownames(sc.dat)

sc.dat.save <- rownames_to_column(as.data.frame(sc.dat), "cell_type")


write_csv(as.data.frame(sc.dat.save), "scRNAseq_from_paper/average_per_cell_from_scRNAseq.csv")

#save count matrix 
gene_dat <- GetAssayData(dat_SUV, assay = "SVA_batchCorrected", layer = "counts")


bk.dat_df <- as.data.frame(bk.dat)
bk.dat_df$gene <- rownames(bk.dat_df)
write_csv(bk.dat_df, "scRNAseq_from_paper/GeoMX_data_counts_SVA_batchCorrected.csv")

# SUV - Create and save a file for each of the options  -----------------------------------------------------------------------------------------------------------------------------

# Create a new directory to store the results
output_dir <- "clean_data_manuscript/instaprism_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

## Loop through each dataset in the list and run InstaPrism
for (i in seq_along(dat_to_run_SUV)) {
  dat_name <- names(dat_to_run_SUV)[i]  # Get the dataset name
  dat_plate <- dat_to_run_SUV[[i]]
  
  # Prep GeoMX dataset for InstaPrism
  gene_dat <- GetAssayData(dat_plate, assay = "SVA_batchCorrected", slot = "counts")
  bk.dat <- as.matrix(t(gene_dat))
  
  # Single-cell data
  sc.dat <- AverageExpression(RDS, group.by = "merged_type") # Average has the beneift of not being biased by the proportion of cells in the original dataset. 
  sc.dat <- t(as.matrix(sc.dat$RNA))
  
  # Labels
  cell.type.labels <- rownames(sc.dat)
  cell.state.labels <- rownames(sc.dat)
  
  # Construct prism object
  myPrism <- new.prism(
    reference = sc.dat,
    mixture = bk.dat,
    input.type = "count.matrix",
    key = NULL,
    cell.type.labels = cell.type.labels,
    cell.state.labels = cell.state.labels
  )
  
  # Run InstaPrism
  InstaPrism.res.initial <- InstaPrism(input_type = 'prism', prismObj = myPrism, n.core = 16)
  
  # Save the InstaPrism object with the dataset name in the file name
  output_file <- file.path(output_dir, paste0("InstaPrism_result_", dat_name, ".rdata"))
  save(InstaPrism.res.initial, file = output_file)
  
  cat("Saved InstaPrism result for dataset", dat_name, "to", output_file, "\n")
}






