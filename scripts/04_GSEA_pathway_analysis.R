#' ==============================================================================
#' Differential Expression to GSEA Analysis
#' ==============================================================================
#' Author: Caroline Duncombe
#' Date: 2024-10-24
#' Description: Generate GSEA results for Hallmark and KEGG gene sets
#' ==============================================================================

# SETUP ------------------------------------------------------------------------

# Load required libraries
library(SEARchways)
library(BIGpicture)
library(dplyr)
library(tibble)
library(readr)
library(here)

# Clear environment
rm(list = ls())

# Define directories
input_dir <- here("clean_data", "DE_results")
output_dir <- here("clean_data", "GSEA")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Analysis parameters
set.seed(123)
N_PERMUTATIONS <- 100000
N_PROCESSORS <- 4

# LOAD DATA --------------------------------------------------------------------

# Load interaction model results
DE_interaction <- readRDS(file.path(input_dir, "DE_interaction_model_all_contrasts.rds")) %>%
  bind_rows(.id = "contrast")

# Load contrast model results - all contrasts
DE_contrast_all <- readRDS(file.path(input_dir, "DE_contrast_model_all_coefficients.rds")) %>%
  bind_rows(.id = "contrast")

# Load contrast model results - selected contrasts
DE_contrast_selected <- readRDS(file.path(input_dir, "DE_contrast_model_specific_contrasts.rds")) %>%
  bind_rows(.id = "contrast_name") %>%
  mutate(
    contrast_name = sub("^Group_Type_", "", contrast_name),
    contrast_name = sub("_differential_expression$", "", contrast_name)
  )

########################################################################################################################################
# STEP 2: RUN GSEA
########################################################################################################################################

# FOR INTERACTION MODEL TERMS 

set.seed(123)

ordered_listall <- DE_interaction %>%  # Include pvalue for use in mutation
  mutate(metric_for_use = log2FoldChange) %>% 
  select(contrast, gene, metric_for_use) %>%# Use log10 for p-values
  arrange(contrast, metric_for_use)  # Arrange by contrast and the new metric

# Question is the symbol a useful fromm here?
#Hallmark
gsea_listall_H <- BIGsea(gene_df=ordered_listall, category="H", ID="SYMBOL", species = "mouse", processors = 4, rand = "multi", nperm = 100000) %>% mutate(database = "hallmark")

gsea_listall_K <- BIGsea(gene_df=ordered_listall, category="C2", subcategory = "CP:KEGG", ID="SYMBOL", species = "mouse", processors = 4, rand = "multi", nperm = 100000) %>% mutate(database = "kegg") 

library(dplyr)  
gsea_listall_K <- gsea_listall_K %>% select(!gs_subcat)

GSEA_output_H_K <- rbind(gsea_listall_H,gsea_listall_K)

saveRDS(GSEA_output_H_K, file = "clean_data/GSEA/GSEA_output_H_K_100000perm_interaction_model.rds")

########################################################################################################################################

# FOR ALL CONTRASTS
set.seed(123)

ordered_listall <- DE_contrast_all %>%  # Include pvalue for use in mutation
  mutate(metric_for_use = log2FoldChange) %>% 
  select(contrast, gene, metric_for_use) %>%# Use log10 for p-values
  arrange(contrast, metric_for_use)  # Arrange by contrast and the new metric


#Hallmark
gsea_listall_H <- BIGsea(gene_df=ordered_listall, category="H", ID="SYMBOL", species = "mouse", processors = 4, rand = "multi", nperm = 100000) %>% mutate(database = "hallmark")

gsea_listall_K <- BIGsea(gene_df=ordered_listall, category="C2", subcategory = "CP:KEGG", ID="SYMBOL", species = "mouse", processors = 4, rand = "multi", nperm = 100000) %>% mutate(database = "kegg") 

library(dplyr)  
gsea_listall_K <- gsea_listall_K %>% select(!gs_subcat)

GSEA_output_H_K <- rbind(gsea_listall_H,gsea_listall_K)

saveRDS(GSEA_output_H_K, file = "clean_data/GSEA/GSEA_output_H_K_100000perm_contrast_all_terms.rds")


########################################################################################################################################

# FOR SELECTED CONTRAST TERMS

set.seed(123)

ordered_listall <- DE_contrast %>%  # Include pvalue for use in mutation
  mutate(metric_for_use = log2FoldChange) %>% 
  select(contrast_name, gene, metric_for_use) %>%# Use log10 for p-values
  arrange(contrast_name, metric_for_use)  # Arrange by contrast and the new metric

#Hallmark
gsea_listall_H <- BIGsea(gene_df=ordered_listall, category="H", ID="SYMBOL", species = "mouse", processors = 4, rand = "multi", nperm = 100000) %>% mutate(database = "hallmark")

gsea_listall_K <- BIGsea(gene_df=ordered_listall, category="C2", subcategory = "CP:KEGG", ID="SYMBOL", species = "mouse", processors = 4, rand = "multi", nperm = 100000) %>% mutate(database = "kegg") 

gsea_listall_K <- gsea_listall_K %>% select(!gs_subcat)

GSEA_output_H_K <- rbind(gsea_listall_H,gsea_listall_K)

table(GSEA_output_H_K$group)

saveRDS(GSEA_output_H_K, file = "clean_data/GSEA/GSEA_output_H_K_100000perm_contrasts_selected_term.rds")

###### ALL DONE #########


