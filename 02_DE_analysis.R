#' ---
#' title: "Differential Expression Analysis with DESeq2"
#' author: "Caroline Duncombe"
#' date: "`r Sys.Date()`"
#' description: "DESeq2 analysis on GeoMX data using SVA batch-corrected counts"
#' ---

# ==============================================================================
# SETUP
# ==============================================================================

# Load required libraries
library(here)
library(GeomxTools)
library(DESeq2)
library(Seurat)
library(dplyr)
library(openxlsx)

# Clear environment
rm(list = ls())

# Set output directory
output_dir <- "clean_data/DE_results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# LOAD DATA
# ==============================================================================

# Load normalized, batch-corrected Seurat object
seuratObj <- readRDS("clean_data/srt/QC_normalized_data_200um_SVA.rds")

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' Create DESeq2 dataset from Seurat object
#'
#' @param seurat_obj Seurat object with SVA_batchCorrected assay
#' @param design_formula Formula for DESeq2 design matrix
#' @param ref_levels Named list of reference levels for factors (e.g., list(Liver_zone = "PP", Group = "M"))
#' @return DESeqDataSet object
#' 
create_dds <- function(seurat_obj, design_formula, ref_levels = NULL) {
  
  # Extract counts and metadata
  cts <- Seurat::GetAssayData(seurat_obj, assay = "SVA_batchCorrected", layer = "counts")
  coldata <- seurat_obj@meta.data
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(cts), 0),
    colData = coldata,
    design = design_formula
  )
  
  # Set reference levels if provided
  if (!is.null(ref_levels)) {
    for (factor_name in names(ref_levels)) {
      dds[[factor_name]] <- relevel(dds[[factor_name]], ref = ref_levels[[factor_name]])
    }
  }
  
  # Run DESeq2
  dds <- DESeq(dds, minReplicatesForReplace = Inf)
  
  return(dds)
}


#' Extract and shrink results for a specific contrast
#'
#' @param dds DESeqDataSet object
#' @param contrast Character vector of length 3 (factor, numerator, denominator) or coefficient name
#' @param shrink_method Method for LFC shrinkage ("apeglm" or "ashr")
#' @return Data frame with differential expression results
#' 
extract_results <- function(dds, contrast, shrink_method = "ashr") {
  
  # Get results based on contrast type
  if (length(contrast) == 1) {
    # Coefficient name provided
    res <- results(dds, name = contrast, cooksCutoff = FALSE, independentFiltering = FALSE)
    res_shrink <- lfcShrink(dds, coef = contrast, type = "apeglm")
  } else {
    # Contrast vector provided
    res <- results(dds, contrast = contrast, cooksCutoff = FALSE, independentFiltering = FALSE)
    res_shrink <- lfcShrink(dds, contrast = contrast, type = shrink_method)
  }
  
  # Convert to data frame
  res_df <- as.data.frame(res_shrink)
  res_df$gene <- rownames(res_df)
  rownames(res_df) <- NULL
  
  return(res_df)
}

#' Save results list to Excel file with each contrast as a separate sheet
#'
#' @param results_list Named list of data frames, one per contrast
#' @param output_path Path for Excel file
#' @param max_sheet_name_length Maximum length for sheet names (default: 31, Excel limit)
#' 
save_results_to_excel <- function(results_list, output_path, max_sheet_name_length = 31) {
  
  # Create workbook
  wb <- createWorkbook()
  
  # Add each contrast as a sheet
  for (contrast_name in names(results_list)) {
    
    # Truncate sheet name if needed (Excel has 31 character limit)
    sheet_name <- contrast_name
    if (nchar(sheet_name) > max_sheet_name_length) {
      sheet_name <- substr(sheet_name, 1, max_sheet_name_length)
    }
    
    # Make sure sheet name is unique
    sheet_name <- make.names(sheet_name, unique = TRUE)
    
    # Add worksheet
    addWorksheet(wb, sheetName = sheet_name)
    
    # Write data
    writeData(wb, sheet = sheet_name, x = results_list[[contrast_name]])
    
    # Auto-size columns
    setColWidths(wb, sheet = sheet_name, cols = 1:ncol(results_list[[contrast_name]]), widths = "auto")
  }
  
  # Save workbook
  saveWorkbook(wb, output_path, overwrite = TRUE)
  message("Saved Excel file: ", output_path)
}

save_results_to_rds <- function(results_list, output_path) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(results_list, output_path)
  message("Saved RDS list: ", output_path)
}


# ==============================================================================
# MODEL 1: INTERACTION MODEL
# ==============================================================================
# Design: Group*Type + Group:Liver_zone + Liver_zone
# Includes all samples (mock, schizont, and bystander)
# ==============================================================================

# Create DESeq2 dataset
dds_interaction <- create_dds(
  seurat_obj = seuratObj,
  design_formula = ~ Group * Type + Group:Liver_zone + Liver_zone,
  ref_levels = list(Liver_zone = "PP", Group = "M")
)

# Extract results for all coefficients
interaction_contrasts <- resultsNames(dds_interaction)
interaction_contrasts <- interaction_contrasts[!grepl("^Intercept", interaction_contrasts)]

# Initialize results list
interaction_results <- list()

for (contrast_name in interaction_contrasts) {
  tryCatch({
    de_res <- extract_results(dds_interaction, contrast = contrast_name, shrink_method = "apeglm")
    de_res$contrast <- contrast_name
    interaction_results[[contrast_name]] <- de_res
  }, error = function(e) {
    warning("Failed to process contrast ", contrast_name, ": ", e$message)
  })
}

# Combine all results
interaction_combined <- bind_rows(interaction_results)

save_results_to_excel(interaction_results,file.path(output_dir, "DE_interaction_model_all_contrasts.xlsx"))

save_results_to_rds(interaction_results,file.path(output_dir, "DE_interaction_model_all_contrasts.rds"))

# Save DESeqDataSet object
saveRDS(dds_interaction, file.path(output_dir, "dds_interaction_model.rds"))


# ==============================================================================
# MODEL 2: CONTRAST MODEL (NO BYSTANDER)
# ==============================================================================
# Design: Group_Type + Liver_zone
# Excludes bystander samples, compares specific group-type combinations
# ==============================================================================

# Filter out bystander samples
seuratObj_noby <- subset(seuratObj, subset = Type != "bystander")

# Create Group_Type combined variable
seuratObj_noby@meta.data <- seuratObj_noby@meta.data %>%
  mutate(Group_Type = paste0(Group, "_", Type))

# Create DESeq2 dataset
cts_noby <- Seurat::GetAssayData(seuratObj_noby, assay = "SVA_batchCorrected", layer = "counts")
coldata_noby <- seuratObj_noby@meta.data

dds_contrast <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(cts_noby), 0),
  colData = coldata_noby,
  design = ~ Group_Type + Liver_zone
)

# Set reference levels
dds_contrast$Liver_zone <- relevel(dds_contrast$Liver_zone, ref = "PP")
dds_contrast$Group_Type <- relevel(dds_contrast$Group_Type, ref = "M_mock")

# Run DESeq2
dds_contrast <- DESeq(dds_contrast, minReplicatesForReplace = Inf)

# ------------------------------------------------------------------------------
# Extract results for all standard coefficients
# ------------------------------------------------------------------------------

contrast_coefficients <- resultsNames(dds_contrast)
contrast_coefficients <- contrast_coefficients[!grepl("^Intercept", contrast_coefficients)]

contrast_results_all <- list()

for (coef_name in contrast_coefficients) {
  tryCatch({
    de_res <- extract_results(dds_contrast, contrast = coef_name, shrink_method = "apeglm")
    de_res$contrast <- coef_name
    contrast_results_all[[coef_name]] <- de_res
  }, error = function(e) {
    warning("Failed to process coefficient ", coef_name, ": ", e$message)
  })
}

# Combine and save all coefficient results
contrast_combined_all <- bind_rows(contrast_results_all)
save_results_to_excel(contrast_results_all,file.path(output_dir, "DE_contrast_model_all_coefficients.xlsx"))
save_results_to_rds(contrast_results_all,file.path(output_dir, "DE_contrast_model_all_coefficients.rds"))


# ------------------------------------------------------------------------------
# Extract specific contrasts of biological interest
# ------------------------------------------------------------------------------

# Define specific contrasts for manuscript
specific_contrasts <- list(
  F_schizont_vs_F_mock = c("Group_Type", "F_schizont", "F_mock"),
  M_schizont_vs_M_mock = c("Group_Type", "M_schizont", "M_mock"),
  ORX_schizont_vs_ORX_mock = c("Group_Type", "ORX_schizont", "ORX_mock"),
  F_mock_vs_M_mock = c("Group_Type", "F_mock", "M_mock"),
  ORX_mock_vs_M_mock = c("Group_Type", "ORX_mock", "M_mock"),
  F_mock_vs_ORX_mock = c("Group_Type", "F_mock", "ORX_mock")
)

# Extract and save each specific contrast
specific_results <- list()

for (contrast_name in names(specific_contrasts)) {
  contrast_vec <- specific_contrasts[[contrast_name]]
  
  tryCatch({
    de_res <- extract_results(dds_contrast, contrast = contrast_vec, shrink_method = "ashr")
    de_res$contrast <- contrast_name
    specific_results[[contrast_name]] <- de_res
    
    # Save individual file
    save_de_results(
      de_res, 
      file.path(output_dir, "specific_contrasts", paste0(contrast_name, "_DE.csv"))
    )
    
  }, error = function(e) {
    warning("Failed to process contrast ", contrast_name, ": ", e$message)
  })
}

# Combine all specific contrasts
specific_combined <- bind_rows(specific_results)

save_results_to_excel(specific_results,file.path(output_dir, "DE_contrast_model_specific_contrasts.xlsx"))
save_results_to_rds(specific_results,file.path(output_dir, "DE_contrast_model_specific_contrasts.rds"))


# Save DESeqDataSet object
saveRDS(dds_contrast, file.path(output_dir, "dds_contrast_model.rds"))

