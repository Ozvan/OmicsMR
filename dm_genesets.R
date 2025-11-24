# Load necessary libraries
library(tidyverse)
library(readxl)
library(biomaRt)
library(data.table)

# ==============================================================================
# 1. PROCESS ALL TESTED GENES
# ==============================================================================

# Helper function to convert nested list structure to a flat, named list of Gene vectors
process_tested_data <- function(data_list, prefix) {
  # Iteratively process each data frame in the list
  purrr::imap(data_list, ~ {
    df <- as.data.frame(.x)
    
    # Check if the data frame has a column (common in single-column RData exports)
    if (ncol(df) == 1) {
      old_col_name <- names(df)[1] 
      
      df %>% 
        dplyr::rename(Gene = !!rlang::sym(old_col_name)) %>% # Rename the single column to 'Gene'
        dplyr::select(Gene) %>%
        dplyr::mutate(dataset = paste0(prefix, .y))
    } else {
      # Handle cases where data frames might have multiple columns, defaulting to a key column name
      # This part assumes that the gene list is stored in the first column if no rename is specified.
      warning(paste("Data frame in", prefix, .y, "has multiple columns. Using first column as 'Gene'."))
      df %>% 
        dplyr::rename(Gene = 1) %>% 
        dplyr::select(Gene) %>%
        dplyr::mutate(dataset = paste0(prefix, .y))
    }
  })
}

# Load the RData files
load("AllTestedGenes_eQTL_pQTL_Blood_ENSG.RData")
load("AllTestedGenes_eQTL_TissueSpecific.RData")

# --- Process Blood eQTL/pQTL ---
blood_eqtl_list <- process_tested_data(Tested.eQTL.pQTL.blood.ENSG$eQTL, "blood.eqtl.")
blood_pqtl_list <- process_tested_data(Tested.eQTL.pQTL.blood.ENSG$pQTL, "blood.pqtl.")

# --- Process Tissue eQTL ---
tissue_names <- c(
  "Subcutaneous Adipose" = "sat", 
  "Visceral Adipose" = "vat", 
  "Skeletal muscle" = "muscle", 
  "Liver" = "liver", 
  "Brain Hypothalamus" = "brain", 
  "Pancreas" = "pancreas", 
  "Pancreatic Islets" = "islets"
)
tissue_eqtl_list <- process_tested_data(res.alltested.alltissues, "tissue.")
names(tissue_eqtl_list) <- tissue_names[names(tissue_eqtl_list)] # Use clearer short names

# Combine all tested lists
tested_data_list <- c(blood_eqtl_list, blood_pqtl_list, tissue_eqtl_list)
all_tested_df <- bind_rows(tested_data_list)
all_tested_ids <- unique(all_tested_df$Gene)

# Clean up environment
rm(Tested.eQTL.pQTL.blood.ENSG, res.alltested.alltissues, blood_eqtl_list, blood_pqtl_list, tissue_eqtl_list)
rm(all_tested_df, tissue_names, tested_data_list)

# ==============================================================================
# 2. PROCESS SIGNIFICANT GENES (Q-value < 0.05)
# ==============================================================================

# Define column names for the large MR excel sheets (moved for clarity)
# NOTE: These vectors must exactly match the number and order of columns in 'media-3.xlsx' sheets.

# Helper function to read, rename, and filter the MR tables
process_significant_mr <- function(filepath, sheet_num, col_names, id_col) {
  df <- read_excel(filepath, sheet_num, skip = 1)
  
  if (ncol(df) != length(col_names)) {
    warning(paste("Column count mismatch in sheet", sheet_num, ". Expected", length(col_names), "got", ncol(df)))
  }
  
  colnames(df) <- col_names
  
  df %>% 
    dplyr::filter(is.na(remove_duplicate)) %>%
    dplyr::mutate(across(matches("_qvalue$"), as.numeric)) %>%
    # Filter for significance in ANY q-value column (non-NA and < 0.05)
    dplyr::filter(if_any(ends_with("_qvalue"), ~ !is.na(.) & . < 0.05)) %>%
    dplyr::pull({{id_col}}) # Dynamically pull the gene/protein ID column
}

# --- Blood eQTL ---
BLOOD_EQTL_COLS <- c(
  "mol.trait", "hgnc_symbol", "metaG_nstudies", "metaG_beta", "metaG_se", "metaG_pvalue", "metaG_qvalue", "metaG_I2",
  "GenEUR_method", "GenEUR_nsnp", "GenEUR_beta", "GenEUR_se", "GenEUR_pvalue", "GenEUR_qvalue", "GenEUR_PASSMRsensitivity", "GenEUR_PWCoCoH4",
  "GenoaAFR_method", "GenoaAFR_nsnp", "GenoaAFR_beta", "GenoaAFR_se", "GenoaAFR_pvalue", "GenoaAFR_qvalue", "GenoaAFR_PASSMRsensitivity", "GenoaAFR_PWCoCoH4",
  "GalaPR_method", "GalaPR_nsnp", "GalaPR_beta", "GalaPR_se", "GalaPR_pvalue", "GalaPR_qvalue", "GalaPR_PASSMRsensitivity", "GalaPR_PWCoCoH4", "remove_duplicate"
)
significant_blood_eqtl_ids <- process_significant_mr("media-3.xlsx", 2, BLOOD_EQTL_COLS, id_col = "mol.trait")
# mol.trait contains the ENSG ID for eQTL

# --- Blood pQTL ---
BLOOD_PQTL_COLS <- c(
  "protein", "hgnc_symbol", "metaP_nstudies", "metaP_beta", "metaP_se", "metaP_pvalue", "metaP_qvalue", "metaP_I2",
  "deCODEEUR_method", "deCODEEUR_nsnp", "deCODEEUR_beta", "deCODEEUR_se", "deCODEEUR_pvalue", "deCODEEUR_qvalue", "deCODEEUR_PASSMRsensitivity", "deCODEEUR_PWCoCoH4",
  "ARICAFR_method", "ARICAFR_nsnp", "ARICAFR_beta", "ARICAFR_se", "ARICAFR_pvalue", "ARICAFR_qvalue", "ARICAFR_PASSMRsensitivity", "ARICAFR_PWCoCoH4",
  "ARICAFR2_method", "ARICAFR2_nsnp", "ARICAFR2_beta", "ARICAFR2_se", "ARICAFR2_pvalue", "ARICAFR2_qvalue", "ARICAFR2_PASSMRsensitivity", "ARICAFR2_PWCoCoH4", "remove_duplicate"
)
significant_blood_pqtl_ids <- process_significant_mr("media-3.xlsx", 3, BLOOD_PQTL_COLS, id_col = "hgnc_symbol")

message("Fetching HGNC to ENSG mapping from Ensembl...")
ensembl <- useMart('ENSEMBL_MART_ENSEMBL', host = 'https://www.ensembl.org')
ensembl <- useDataset('hsapiens_gene_ensembl', ensembl)

ense_hgnc_map <- getBM(
  attributes = c('hgnc_symbol', 'ensembl_gene_id'),
  mart = ensembl,
  uniqueRows = TRUE,
  filters='hgnc_symbol', values=significant_blood_pqtl_ids
)
significant_blood_pqtl_ids <- ense_hgnc_map$ensembl_gene_id

# --- Tissue eQTL ---
TISSUE_EQTL_COLS <- c(
  "mol.trait", "hgnc_symbol",
  "SAT_method", "SAT_nsnp", "SAT_beta", "SAT_se", "SAT_pvalue", "SAT_qvalue", "SAT_PASSMRsensitivity", "SAT_PWCoCoH4",
  "VAT_method", "VAT_nsnp", "VAT_beta", "VAT_se", "VAT_pvalue", "VAT_qvalue", "VAT_PASSMRsensitivity", "VAT_PWCoCoH4",
  "Brain_method", "Brain_nsnp", "Brain_beta", "Brain_se", "Brain_pvalue", "Brain_qvalue", "Brain_PASSMRsensitivity", "Brain_PWCoCoH4",
  "Liver_method", "Liver_nsnp", "Liver_beta", "Liver_se", "Liver_pvalue", "Liver_qvalue", "Liver_PASSMRsensitivity", "Liver_PWCoCoH4",
  "Muscle_method", "Muscle_nsnp", "Muscle_beta", "Muscle_se", "Muscle_pvalue", "Muscle_qvalue", "Muscle_PASSMRsensitivity", "Muscle_PWCoCoH4",
  "Pancreas_method", "Pancreas_nsnp", "Pancreas_beta", "Pancreas_se", "Pancreas_pvalue", "Pancreas_qvalue", "Pancreas_PASSMRsensitivity", "Pancreas_PWCoCoH4",
  "Islets_method", "Islets_nsnp", "Islets_beta", "Islets_se", "Islets_pvalue", "Islets_qvalue", "Islets_PASSMRsensitivity", "Islets_PWCoCoH4", "remove_duplicate"
)
tissue_mr_df <- read_excel("/Users/ahuerta/Dropbox/Mis_otros_proyectos/Broad/MR_omics/paper/media-3.xlsx", 7, skip = 1)
colnames(tissue_mr_df) <- TISSUE_EQTL_COLS

# Function to filter tissue data for a single q-value column
filter_tissue_sig <- function(tissue_df, qvalue_col) {
  tissue_df %>% 
    dplyr::filter(is.na(remove_duplicate)) %>%
    dplyr::mutate(across(matches("_qvalue$"), as.numeric)) %>%
    dplyr::filter(!!rlang::sym(qvalue_col) < 0.05) %>%
    dplyr::pull(mol.trait) # Pulls the ENSG ID
}

# Apply filtering to each tissue-specific q-value column
qvalue_cols <- c("SAT_qvalue", "VAT_qvalue", "Brain_qvalue", "Liver_qvalue", "Muscle_qvalue", "Pancreas_qvalue", "Islets_qvalue")
tissue_sig_lists <- purrr::map(qvalue_cols, ~ filter_tissue_sig(tissue_mr_df, .x))
names(tissue_sig_lists) <- gsub("_qvalue", "", qvalue_cols) # Name lists clearly

# Combine all significant ENSG ID lists
significant_ensg_lists <- c(
  list(bloodeqtl = significant_blood_eqtl_ids ,bloodpqtl = significant_blood_pqtl_ids,
  tissueeqtl = tissue_sig_lists
))

significant_tested_ids <- unique(unlist(significant_ensg_lists))

# ==============================================================================
# 3. PROCESS CAUSAL GENES
# ==============================================================================

load("AllSigGenes_eQTL_pQTL_rep_ENSG.RData") # all.significant.genes.eQTL.pQTL.ENSG
load("AllSigGenes_eQTL_TissueSpecific.RData") # res.sig.alltissues

causal_ensg_lists <- list(
  blood.eqtl.eur = all.significant.genes.eQTL.pQTL.ENSG$eQTL$EUR,
  blood.eqtl.afr = all.significant.genes.eQTL.pQTL.ENSG$eQTL$AFR,
  blood.eqtl.amr = all.significant.genes.eQTL.pQTL.ENSG$eQTL$AMR,
  blood.eqtl.meta = all.significant.genes.eQTL.pQTL.ENSG$eQTL$MetaAnalysis,
  blood.pqtl.eur = all.significant.genes.eQTL.pQTL.ENSG$pQTL$EUR,
  blood.pqtl.afr = all.significant.genes.eQTL.pQTL.ENSG$pQTL$AFR,
  blood.pqtl.eas = all.significant.genes.eQTL.pQTL.ENSG$pQTL$EAS,
  blood.pqtl.meta = all.significant.genes.eQTL.pQTL.ENSG$pQTL$MetaAnalysis,
  sat = res.sig.alltissues$`Subcutaneous Adipose`,
  vat = res.sig.alltissues$`Visceral Adipose`,
  brain = res.sig.alltissues$`Brain Hypothalamus`,
  liver = res.sig.alltissues$Liver,
  muscle = res.sig.alltissues$`Skeletal muscle`,
  pancreas = res.sig.alltissues$Pancreas,
  islets = res.sig.alltissues$`Pancreatic Islets`
)
causal_tested_ids <- unique(unlist(causal_ensg_lists))

rm(all.significant.genes.eQTL.pQTL.ENSG, res.sig.alltissues)

# ==============================================================================
# 4. PROCESS DM GENE SETS AND MAP IDs
# ==============================================================================

dm_genes_df <- read_excel("DM_genes.xlsx", 2)
# Create initial HGNC-based gene lists
gene_sets_hgnc <- list(
  HPO = unique(na.omit(unlist(dm_genes_df[, grepl("hpo", names(dm_genes_df))]))),
  Monogenic = unique(na.omit(unlist(dm_genes_df[, grepl("monogenic", names(dm_genes_df))]))),
  Congenital = unique(na.omit(unlist(dm_genes_df[, grepl("congenital", names(dm_genes_df))]))),
  Neonatal = unique(na.omit(unlist(dm_genes_df[, grepl("Neonatal.dm", names(dm_genes_df))]))),
  MousePhenotype = unique(na.omit(unlist(dm_genes_df[, grepl("MP", names(dm_genes_df))]))),
  Mahajan = unique(na.omit(unlist(dm_genes_df[, grepl("mahajan", names(dm_genes_df))]))),
  HuGE = unique(na.omit(unlist(dm_genes_df[, grepl("huge", names(dm_genes_df))]))),
  Mahajan_HuGE = unique(na.omit(unlist(dm_genes_df[, c("mahajan.score5.t2d.gwas", "mahajan.score4.t2d.gwas", "mahajan.score3.t2d.gwas", "huge.score")])))
)
original_hgnc_set_sizes <- sapply(gene_sets_hgnc, length)

# Correction for a specific gene name (done before mapping to ENSG)
gene_sets_hgnc$Congenital <- gene_sets_hgnc$Congenital[gene_sets_hgnc$Congenital != "NKX&-1"]

# --- Map HGNC symbols to Ensembl Gene IDs (ENSG) ---
all_hgnc_symbols <- unique(unlist(gene_sets_hgnc))

# Fetch mapping from Ensembl using biomaRt (may take time)
message("Fetching HGNC to ENSG mapping from Ensembl...")
ensembl <- useMart('ENSEMBL_MART_ENSEMBL', host = 'https://www.ensembl.org')
ensembl <- useDataset('hsapiens_gene_ensembl', ensembl)

ense_hgnc_map_df <- getBM(
  attributes = c('hgnc_symbol', 'ensembl_gene_id'),
  mart = ensembl,
  uniqueRows = TRUE,
  filters='hgnc_symbol', values=all_hgnc_symbols
) %>% 
  dplyr::filter(ensembl_gene_id != "")

dmgene_lists_ensg <- lapply(gene_sets_hgnc, function(hgnc_list) {
  current_hgnc_df <- tibble(hgnc_symbol = hgnc_list)
  
  mapped_df <- current_hgnc_df %>%
    dplyr::left_join(ense_hgnc_map_df, by = "hgnc_symbol")
  
  mapped_df %>%
    dplyr::pull(ensembl_gene_id) %>%
    unique() %>%
    na.omit()
})

# The set of tested genes that were *never* significant
never_significant_ids <- setdiff(all_tested_ids, significant_tested_ids)

# ==============================================================================
# 5. GENERATE OVERLAP COUNTS (IRRESPECTIVE OF TISSUE EXPRESSION)
# ==============================================================================

# Function to count overlaps
count_overlap <- function(target, reference) {
  length(intersect(target, reference))
}

# Apply across all DM gene lists (now in ENSG format)
overlap_counts_unfiltered <- purrr::map_dfr(dmgene_lists_ensg, ~ {
  tibble(
    in_all_tested = count_overlap(.x, all_tested_ids),
    in_any_significant = count_overlap(.x, significant_tested_ids),
    in_causal = count_overlap(.x, causal_tested_ids),
    in_never_significant = count_overlap(.x, never_significant_ids)
  )
}, .id = "Gene_Set")

original_sizes_df <- original_hgnc_set_sizes %>%
  tibble::enframe(name = "Gene_Set", value = "original_hgnc_size")
overlap_counts_unfiltered <- overlap_counts_unfiltered %>%
  dplyr::left_join(original_sizes_df, by = "Gene_Set")

overlap_counts_unfiltered %>%
  mutate(perc_tested=(in_all_tested/original_hgnc_size)*100,
         perc_any_significant_among_tested=(in_any_significant/in_all_tested)*100,
         perc_causal_among_tested=(in_causal/in_all_tested)*100,
         perc_never_significant_among_tested=(in_never_significant/in_all_tested)*100) %>%
  dplyr::select(Gene_Set, perc_tested,perc_never_significant_among_tested,perc_any_significant_among_tested,perc_causal_among_tested)
         
         



