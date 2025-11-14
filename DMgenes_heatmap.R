library(data.table)
library(ggplot2)
library(dplyr)
library(purrr)
library(ComplexHeatmap)
library(circlize)

mahaj_ense <- fread("mahajan+HuGE30_ensemblid_mapping.txt")
mahaj_tissue <- fread("MahajanHuge_expressed1+GTEx_ENSGmappings_27.08.25.txt")
map_ref <- mahaj_tissue[!is.na(`HGNC symbol`) & significant_in_any=="yes",
                        .(`HGNC symbol`,Gene)]

mahaj_tissue[!is.na(`HGNC symbol`) & significant_in_any=="yes",.N]
# [1] 141

map_ref_tis <- mahaj_tissue[!is.na(`HGNC symbol`) & 
                              (blood.eqtl.eur=="significant" | 
                                 sat=="significant" | 
                                 vat=="significant" | 
                                 muscle=="significant" | 
                                 liver=="significant" | 
                                 brain=="significant" | 
                                 pancreas=="significant" | 
                                 islets=="significant"),
                            .(`HGNC symbol`,Gene)]

### all qtl
blood.eqtl.afr <- fread("sensitivity_mr_results_eqtl_GENOA_plasma_AA.txt") |>
  dplyr::mutate(eqtl.afr.z = beta/standard_error,
                mol.trait = sub(".*__", "", sub("\\..*", "", exposure))) |>
  dplyr::select(mol.trait,eqtl.afr.z)

blood.eqtl.amr <- fread("mr_results_eqtl_GALAII_plasma_MX.txt") |>
  dplyr::mutate(eqtl.amr.z = beta/standard_error,
                mol.trait = sub(".*__", "", sub("\\..*", "", exposure))) |>
  dplyr::select(mol.trait,eqtl.amr.z)

blood.eqtl.eur <- fread("original_sensitivity_mr_results_eqtl_eQTLGen_Plasma_EUR.txt") |>
  dplyr::mutate(eqtl.eur.z = beta/standard_error) |>
  dplyr::select(mol.trait,eqtl.eur.z)

blood.eqtl.z <- Reduce(function(x, y) merge(x, y, by = "mol.trait", all = T), 
                       list(blood.eqtl.afr,blood.eqtl.amr,blood.eqtl.eur)) |>
  dplyr::filter(mol.trait %in% map_ref$Gene)

### pQTL
blood.pqtl <- fread("pQTL_meta_analysis_deCODEEUR_ARICAFR_KyotoEAS_randommetapval_Aug272024.tsv")

blood.pqtl.z <- blood.pqtl |>
  dplyr::mutate(pqtl.afr.z = b.ARIC_AFR/se.ARIC_AFR,
                pqtl.eas.z = b.Kyoto_EAS/se.Kyoto_EAS,
                pqtl.eur.z = b.deCODE_EUR/se.deCODE_EUR) |>
  dplyr::select(shortname,pqtl.afr.z,pqtl.eas.z,pqtl.eur.z)

blood.pqtl.ense <- merge(x = map_ref,
                         y = blood.pqtl.z,
                         by.x = "HGNC symbol",
                         by.y = "shortname",
                         all = F) |>
  dplyr::select(Gene,pqtl.afr.z,pqtl.eas.z,pqtl.eur.z)

### Blood qtl overall
blood.qtl <- merge(
  x = blood.eqtl.z,
  y = blood.pqtl.ense,
  by.x = "mol.trait",
  by.y = "Gene",
  all = TRUE
) %>%
  rowwise() %>%
  mutate(
    best_z = {
      z_values <- c(eqtl.afr.z, eqtl.amr.z, eqtl.eur.z,
                    pqtl.afr.z, pqtl.eas.z, pqtl.eur.z)
      if (all(is.na(z_values))) {
        NA_real_
      } else {
        z_values[which.max(abs(z_values))]
      }
    }
  ) %>%
  ungroup() %>%
  dplyr::select(mol.trait, best_z) %>%
  rename(Blood = best_z) %>%
  group_by(mol.trait) %>%
  summarise(Blood = Blood[which.max(abs(Blood))], .groups = "drop")


### Other 7 tissues
files <- list.files(path = "tissue/", pattern = "^sensitivity_.*_EUR\\.txt$")

all.tested.z   <- list()

for (f in files) {
  tissue <- sub("sensitivity_mr_results_(eqtl|eQTL)_(eQTLGen|GTEx|TIGER)_", "",
                sub("_EUR\\.txt$", "", f))
  tissue <- gsub("([a-z])([A-Z])", "\\1 \\2", tissue)
  
  res <- fread(paste0("tissue/",f))
  
  if (tissue %in% c("Pancreatic Islets", "Plasma")) {
    res[,sig := ifelse(p_value_fdr < 0.05 & 
                         !DiffDirection & 
                         !FlagSensitivity & 
                         PWCoCo_unconditioned_H4 > 0.8, 1, 0)]
  } else {
    res[,sig := ifelse(p_value_fdr < 0.05 & 
                         !DiffDirection & 
                         !FlagSensitivity & 
                         PWCoCo_unconditioned_H4 > 0.8, 1, 0)]
  }
  
  res <- res[,.(mol.trait,beta,standard_error,sig)]

  res[, mol.trait := sapply(mol.trait, function(z) strsplit(z, split = "\\.")[[1]][1])]
  res[, z := beta / standard_error]
  
  res.z <- res[mol.trait %in% map_ref$Gene,.(mol.trait, z)]
  
  setnames(res.z, "z", tissue)
  
  all.tested.z[[tissue]]   <- res.z
}

all.tested.z[["Blood"]] <- blood.qtl

all.tissue.z <- Reduce(function(x, y) merge(x, y, by = "mol.trait", all = T), all.tested.z) |>
  unique()

plot.z.hgnc <- merge(x = all.tissue.z,
                     y = map_ref,
                     by.x = "mol.trait",
                     by.y = "Gene",
                     all = F)

plot.z.hgnc <- plot.z.hgnc[order(`HGNC symbol`),.(`HGNC symbol`,
                                                  Blood,
                                                  `Brain Hypothalamus`,
                                                  Liver,
                                                  Pancreas,
                                                  `Pancreatic Islets`,
                                                  `Skeletal Muscle`,
                                                  `Subcutaneous Adipose`,
                                                  `Visceral Adipose`)]
colnames(plot.z.hgnc)[1] <- "Gene/Protein"


## ComplexHeatmap ----
z_mat <- as.matrix(as.data.frame(plot.z.hgnc[, -1]))
rownames(z_mat) <- plot.z.hgnc$`Gene/Protein`

all.test <- fread("alltested_TPM0.1_GTExv10_20.08.25.txt")

sig.indi <- all.test |>
  mutate(blood=ifelse(blood.eqtl.eur=="tested_significant"|
                        blood.eqtl.afr=="tested_significant"|
                        blood.eqtl.amr=="tested_significant"|
                        blood.pqtl.eur=="tested_significant"|
                        blood.pqtl.eas=="tested_significant", "*","")) |>
  dplyr::select(hgnc_symbol,
                blood,
                brain,
                liver,
                pancreas,
                islets,
                muscle,
                sat,
                vat)

sig_df <- merge(x = plot.z.hgnc[,.(`Gene/Protein`)],
                y = sig.indi,
                by.x = "Gene/Protein",
                by.y = "hgnc_symbol",
                all = F)

colnames(sig_df) <- colnames(plot.z.hgnc)
sig_mat <- as.matrix(as.data.frame(sig_df[,-1]))
rownames(sig_mat) <- sig_df$`Gene/Protein`
sig_mat <- cbind(sig_mat[,1],
                 apply(sig_mat[,-1], 2, function(x) ifelse(x == "tested_significant", "*", "")))

col_fun <- colorRamp2(c(-20, 0, 20), c("#2b6cb0", "white", "#c53030"))

#--- main heatmap ---
ht.1 <- Heatmap(
  t(z_mat[1:47,]),  
  name = "Z",
  col = col_fun,
  cluster_rows = F,      
  cluster_columns = F,   
  show_row_names = T,
  show_column_names = T,
  row_names_side = "left",
  column_names_side = "bottom",
  na_col = "grey90",
  column_names_gp = gpar(fontsize = 8,
                         fontface = "bold.italic"),
  border = T, 
  cell_fun = function(j, i, x, y, w, h, fill) {
    lab <- t(sig_mat[1:47,])[i, j]
    if (!is.na(lab)) {
      grid.text(lab, x, y, gp = gpar(fontsize = 12))
    }
  },
  width = unit(17, "in"), 
  height = unit(3, "in")  
)

draw(ht.1)

ht.2 <- Heatmap(
  t(z_mat[48:94,]),  
  name = "Z",
  col = col_fun,
  cluster_rows = F,      
  cluster_columns = F,   
  show_row_names = T,
  show_column_names = T, 
  row_names_side = "left",
  column_names_side = "bottom",
  na_col = "grey90",
  column_names_gp = gpar(fontsize = 8,
                         fontface = "bold.italic"),
  border = T, 
  cell_fun = function(j, i, x, y, w, h, fill) {
    lab <- t(sig_mat[48:94,])[i, j]
    if (!is.na(lab)) {
      grid.text(lab, x, y, gp = gpar(fontsize = 12))
    }
  },
  width = unit(17, "in"),   
  height = unit(3, "in")   
)

draw(ht.2)

ht.3 <- Heatmap(
  t(z_mat[95:141,]),  
  name = "Z",
  col = col_fun,
  cluster_rows = F,     
  cluster_columns = F,  
  show_row_names = T,
  show_column_names = T,
  row_names_side = "left",
  column_names_side = "bottom",
  na_col = "grey90",
  column_names_gp = gpar(fontsize = 8,
                         fontface = "bold.italic"),
  border = T, 
  cell_fun = function(j, i, x, y, w, h, fill) {
    lab <- t(sig_mat[95:141,])[i, j]
    if (!is.na(lab)) {
      grid.text(lab, x, y, gp = gpar(fontsize = 12))
    }
  },
  width = unit(17, "in"),   
  height = unit(3, "in")   
)

draw(ht.3)
