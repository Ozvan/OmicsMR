library(data.table)
library(ggplot2)
library(dplyr)
`%notin%` <- Negate(`%in%`) # negation of the given function: %in%
path="../Comparisons_clusters/Data/"

load(paste0(path, "AllSigGenes_eQTL_TissueSpecific.RData"))
load(paste0(path, "AllSigGenes_eQTL_TissueSpecific_Bysign.RData"))
load(paste0(path, "AllSigGenes_eQTL_pQTL_rep.RData"))
load(paste0(path, "AllTestedGenes_eQTL_TissueSpecific.RData"))
load(paste0(path, "AllTestedGenes_eQTL_pQTL_Blood.RData"))
res.sig.alltissues$Blood <- unlist(all.significant.genes.eQTL.pQTL$eQTL$EUR)
res.alltested.alltissues$Blood <- unlist(Tested.eQTL.pQTL.blood$eQTL$eQTLGen_EUR)


cluster.Suzuki <- read.table(paste0(path, "ST6_Suzuki24_IndexVariants_Clusters.txt"), header = T, as.is = T, sep = "\t")
rownames(cluster.Suzuki) <- cluster.Suzuki$Index.SNV
FUMA.snps <- read.table(paste0(path, "FUMA_job519881/snps.txt"), header = T)
FUMA.snps.genes <- lapply(FUMA.snps$nearestGene, function(z) strsplit(z, split = ":"))
FUMA.snps.genes.df <- data.frame(SNP = rep(FUMA.snps$rsID, unlist(lapply(FUMA.snps.genes, function(z) length(z[[1]])))),
                                 Gene = unlist(FUMA.snps.genes))
FUMA.snps.genes.df$cluster <- cluster.Suzuki[FUMA.snps.genes.df$SNP, "Cluster.assignment"]
Gene.bycluster <- sapply(unique(FUMA.snps.genes.df$cluster), function(z) subset(FUMA.snps.genes.df, cluster == z)$Gene)
cluster_map <- stack(Gene.bycluster)
colnames(cluster_map) <- c("Gene","Cluster")


# hgnc --------------------------------------------------------------------
ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://grch37.ensembl.org")
res.sig.alltissues.hgnc <- res.sig.alltissues

for(tissue in names(res.sig.alltissues)){
  res.sig.alltissues.hgnc[[tissue]] <- biomaRt::getBM(filters= "ensembl_gene_id",
                                                      attributes= c("ensembl_gene_id","hgnc_symbol"),
                                                      values=res.sig.alltissues[[tissue]],
                                                      mart= ensembl)$hgnc_symbol

  res.sig.alltissues.hgnc[[tissue]] <- res.sig.alltissues.hgnc[[tissue]][which(res.sig.alltissues.hgnc[[tissue]] != "")]
}

res.alltested.alltissues.hgnc <- res.alltested.alltissues
for(tissue in names(res.alltested.alltissues.hgnc)){
  res.alltested.alltissues.hgnc[[tissue]] <- biomaRt::getBM(filters= "ensembl_gene_id", 
                                                            attributes= c("ensembl_gene_id","hgnc_symbol"),
                                                            values=res.alltested.alltissues[[tissue]],
                                                            mart= ensembl)$hgnc_symbol

  res.alltested.alltissues.hgnc[[tissue]] <- res.alltested.alltissues.hgnc[[tissue]][which(res.alltested.alltissues.hgnc[[tissue]] != "")]
}

# res.sig.alltissues.hgnc <- lapply(names(res.sig.alltissues.bysign.hgnc$beta_neg), function(z) unique(c(res.sig.alltissues.bysign.hgnc$beta_neg[[z]], res.sig.alltissues.bysign.hgnc$beta_pos[[z]])))
# names(res.sig.alltissues.hgnc) <- names(res.sig.alltissues.bysign.hgnc$beta_neg)

# All sig -----------------------------------------------------------------
results <- list()

for (tissue in names(res.sig.alltissues.hgnc)) {
  sig_genes <- res.sig.alltissues.hgnc[[tissue]]
  tested_genes <- res.alltested.alltissues.hgnc[[tissue]]
  
  bg_genes <- intersect(tested_genes, unique(cluster_map$Gene))
  if (length(bg_genes) == 0) next  
  
  sig_genes <- intersect(sig_genes, bg_genes)
  
  for (cluster in unique(cluster_map$Cluster)) {
    cluster_genes <- intersect(cluster_map$Gene[cluster_map$Cluster == cluster], bg_genes)
    if (length(cluster_genes) == 0) next  
    
    a <- sum(sig_genes %in% cluster_genes)
    b <- sum(cluster_genes %notin% sig_genes)
    c <- sum(sig_genes %notin% cluster_genes)
    d <- length(setdiff(bg_genes, union(sig_genes, cluster_genes)))
    
    mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    dimnames(mat) <- list(Sig=c("Y","N"),
                          Cluster=c("Y","N"))
    
    if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) next
    
    test <- fisher.test(mat, simulate.p.value = TRUE)
    
    results[[paste(tissue, cluster, sep = "_")]] <- data.frame(
      Tissue = tissue,
      Cluster = cluster,
      Overlap = a,
      OddsRatio = test$estimate,
      Pvalue = test$p.value
    )
  }
}


res_df <- do.call(rbind, results)
res_df_no0 <- res_df %>% filter(Overlap > 0)
length(which(res_df_no0$Pvalue < 0.05))
# [1] 4

## FDR (across all)
res_df_no0$FDR <- p.adjust(res_df_no0$Pvalue, method = "fdr")
length(which(res_df_no0$FDR < 0.05))
# [1] 0  No significant results after FDR across all p-values

res_df_no0$negLog10P <- -log10(res_df_no0$Pvalue)

x_limits <- range(res_df_no0$OddsRatio, na.rm = TRUE)

ggplot(res_df_no0, aes(x = OddsRatio, 
                         y = negLog10P, 
                         color = Cluster, 
                         shape = Tissue)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    x = "Odds Ratio",
    y = "-log10(P)",
    # title = "Based on all significant genes"
  ) +
  scale_x_continuous(
    limits = x_limits,
    sec.axis = sec_axis(
      trans = ~., 
      breaks = c(x_limits[1], x_limits[2]),
      labels = c("Depleted", "Enriched")
    )
  ) +
  theme_minimal() +
  theme(
    axis.text.x.top = element_text(face = "bold", color = "black"),
    axis.ticks.x.top = element_blank(),
    axis.title.x.top = element_blank()
  ) +
  scale_shape_manual(values = c(10,15,16,12,9,17,14,8,11))


# By sign -----------------------------------------------------------------
res.sig.alltissues.bysign.hgnc <- res.sig.alltissues.bysign

for(tissue in names(res.sig.alltissues.bysign$beta_neg)){
  res.sig.alltissues.bysign.hgnc$beta_neg[[tissue]] <- biomaRt::getBM(filters= "ensembl_gene_id", 
                                                                      attributes= c("ensembl_gene_id","hgnc_symbol"),
                                                                      values=res.sig.alltissues.bysign$beta_neg[[tissue]],
                                                                      mart= ensembl)$hgnc_symbol

  res.sig.alltissues.bysign.hgnc$beta_neg[[tissue]] <- res.sig.alltissues.bysign.hgnc$beta_neg[[tissue]][which(res.sig.alltissues.bysign.hgnc$beta_neg[[tissue]] != "")]

  res.sig.alltissues.bysign.hgnc$beta_pos[[tissue]] <- biomaRt::getBM(filters= "ensembl_gene_id", 
                                                                      attributes= c("ensembl_gene_id","hgnc_symbol"),
                                                                      values=res.sig.alltissues.bysign$beta_pos[[tissue]],
                                                                      mart= ensembl)$hgnc_symbol

  res.sig.alltissues.bysign.hgnc$beta_pos[[tissue]] <- res.sig.alltissues.bysign.hgnc$beta_pos[[tissue]][which(res.sig.alltissues.bysign.hgnc$beta_pos[[tissue]] != "")]
}

res.anc.EUR.eQTLGen <- fread("sensitivity_mr_results_eqtl_eQTLGen_Blood_EUR.txt")
res.anc.EUR.eQTLGen[p_value_ivw_fdr < 0.05 &
                      !DiffDirection & 
                      !FlagSensitivity & 
                      PWCoCo_unconditioned_H4>0.8, .N]

res.anc.EUR.eQTLGen[, Sig := ifelse(p_value_ivw_fdr < 0.05 &
                                      !DiffDirection &
                                      !FlagSensitivity &
                                      PWCoCo_unconditioned_H4>0.8,
                                    TRUE, FALSE)]
res.anc.EUR.eQTLGen.sig <- res.anc.EUR.eQTLGen[Sig==T]

eQTLGen.sig.hgnc <- biomaRt::getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = res.anc.EUR.eQTLGen.sig$mol.trait,
  mart = ensembl
)$hgnc_symbol

eQTLGen.sig.hgnc.uniq <- unique(eQTLGen.sig.hgnc[eQTLGen.sig.hgnc != "" & !is.na(eQTLGen.sig.hgnc)])

eQTLGen.pos.hgnc <- biomaRt::getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = res.anc.EUR.eQTLGen.sig[beta>0,mol.trait],
  mart = ensembl
)$hgnc_symbol

eQTLGen.pos.hgnc.uniq <- unique(eQTLGen.pos.hgnc[eQTLGen.pos.hgnc != "" & !is.na(eQTLGen.pos.hgnc)])

eQTLGen.neg.hgnc <- biomaRt::getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = res.anc.EUR.eQTLGen.sig[beta<0,mol.trait],
  mart = ensembl
)$hgnc_symbol

eQTLGen.neg.hgnc.uniq <- unique(eQTLGen.neg.hgnc[eQTLGen.neg.hgnc != "" & !is.na(eQTLGen.neg.hgnc)])

res.sig.alltissues.bysign.hgnc$beta_pos[["Blood"]] <- eQTLGen.pos.hgnc.uniq
res.sig.alltissues.bysign.hgnc$beta_neg[["Blood"]] <- eQTLGen.neg.hgnc.uniq

res.alltested.alltissues.hgnc <- res.alltested.alltissues
for(tissue in names(res.alltested.alltissues.hgnc)){
  res.alltested.alltissues.hgnc[[tissue]] <- biomaRt::getBM(filters= "ensembl_gene_id", 
                                                            attributes= c("ensembl_gene_id","hgnc_symbol"),
                                                            values=res.alltested.alltissues[[tissue]],
                                                            mart= ensembl)$hgnc_symbol
  
  res.alltested.alltissues.hgnc[[tissue]] <- res.alltested.alltissues.hgnc[[tissue]][which(res.alltested.alltissues.hgnc[[tissue]] != "")]
}

res.sig.alltissues.hgnc <- lapply(names(res.sig.alltissues.bysign.hgnc$beta_neg), function(z) unique(c(res.sig.alltissues.bysign.hgnc$beta_neg[[z]], res.sig.alltissues.bysign.hgnc$beta_pos[[z]])))
names(res.sig.alltissues.hgnc) <- names(res.sig.alltissues.bysign.hgnc$beta_neg)

# Sig postive -----------------------------------------------------------------
results_pos <- list()

for (tissue in names(res.sig.alltissues.hgnc)) {
  pos_genes <- res.sig.alltissues.bysign.hgnc$beta_pos[[tissue]]
  tested_genes <- res.alltested.alltissues.hgnc[[tissue]]
  
  bg_genes <- intersect(tested_genes, unique(cluster_map$Gene))
  if (length(bg_genes) == 0) next  
  
  pos_genes <- intersect(pos_genes, bg_genes)
  
  for (cluster in unique(cluster_map$Cluster)) {
    cluster_genes <- intersect(cluster_map$Gene[cluster_map$Cluster == cluster], bg_genes)
    if (length(cluster_genes) == 0) next  

    a <- sum(pos_genes %in% cluster_genes)
    b <- sum(cluster_genes %notin% pos_genes)
    c <- sum(pos_genes %notin% cluster_genes)
    d <- length(setdiff(bg_genes, union(pos_genes, cluster_genes)))
    
    mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    dimnames(mat) <- list(Sig=c("Y","N"),
                          Cluster=c("Y","N"))
    
    if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) next
    
    test <- fisher.test(mat, simulate.p.value = TRUE)
    
    results_pos[[paste(tissue, cluster, sep = "_")]] <- data.frame(
      Tissue = tissue,
      Cluster = cluster,
      Overlap = a,
      OddsRatio = test$estimate,
      Pvalue = test$p.value
    )
  }
}

res_df_pos <- do.call(rbind, results_pos)
res_df_pos_no0 <- res_df_pos %>% filter(Overlap > 0)
length(which(res_df_pos_no0$Pvalue < 0.05))
# [1] 2

res_df_pos_no0$FDR <- p.adjust(res_df_pos_no0$Pvalue, method = "fdr")
length(which(res_df_pos_no0$FDR < 0.05))
# [1] 0  No significant results after FDR across all p-values

res_df_pos_no0$negLog10P <- -log10(res_df_pos_no0$Pvalue)

x_limits <- range(res_df_pos_no0$OddsRatio, na.rm = TRUE)

ggplot(res_df_pos_no0, aes(x = OddsRatio, y = negLog10P, color = Cluster, shape = Tissue)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    x = "Odds Ratio",
    y = "-log10(P)",
    # title = "Based on significant genes with a positive effect"
  ) +
  scale_x_continuous(
    limits = x_limits,
    sec.axis = sec_axis(
      trans = ~., 
      breaks = c(x_limits[1], x_limits[2]),
      labels = c("Depleted", "Enriched")
    )
  ) +
  theme_minimal() +
  theme(
    axis.text.x.top = element_text(face = "bold", color = "black"),
    axis.ticks.x.top = element_blank(),
    axis.title.x.top = element_blank()
  ) +
  scale_shape_manual(values = c(10,15,16,12,9,17,14,8,11))


# Sig neg -----------------------------------------------------------------
results_neg <- list()

for (tissue in names(res.sig.alltissues.hgnc)[1:8]) {
  neg_genes <- res.sig.alltissues.bysign.hgnc$beta_neg[[tissue]]
  tested_genes <- res.alltested.alltissues.hgnc[[tissue]]
  
  bg_genes <- intersect(tested_genes, unique(cluster_map$Gene))
  if (length(bg_genes) == 0) next 
  
  neg_genes <- intersect(neg_genes, bg_genes)
  
  for (cluster in unique(cluster_map$Cluster)) {
    cluster_genes <- intersect(cluster_map$Gene[cluster_map$Cluster == cluster], bg_genes)
    if (length(cluster_genes) == 0) next  
    
    a <- sum(neg_genes %in% cluster_genes)
    b <- sum(cluster_genes %notin% neg_genes)
    c <- sum(neg_genes %notin% cluster_genes)
    d <- length(setdiff(bg_genes, union(neg_genes, cluster_genes)))
    
    mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    dimnames(mat) <- list(Sig=c("Y","N"),
                          Cluster=c("Y","N"))
    
    if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) next
    
    test <- fisher.test(mat, simulate.p.value = TRUE)
    
    results_neg[[paste(tissue, cluster, sep = "_")]] <- data.frame(
      Tissue = tissue,
      Cluster = cluster,
      Overlap = a,
      OddsRatio = test$estimate,
      Pvalue = test$p.value
    )
  }
}

res_df_neg <- do.call(rbind, results_neg)
res_df_neg_no0 <- res_df_neg %>% filter(Overlap > 0)
length(which(res_df_neg_no0$Pvalue < 0.05))
# [1] 1

res_df_neg_no0$FDR <- p.adjust(res_df_neg_no0$Pvalue, method = "fdr")
length(which(res_df_neg_no0$FDR < 0.05))
# [1] 0  No significant results after FDR across all p-values

res_df_neg_no0$negLog10P <- -log10(res_df_neg_no0$Pvalue)

x_limits <- range(res_df_neg_no0$OddsRatio, na.rm = TRUE)

ggplot(res_df_neg_no0, aes(x = OddsRatio, y = negLog10P, color = Cluster, shape = Tissue)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    x = "Odds Ratio",
    y = "-log10(P)"
    # ,title = "Based on significant genes with a negative effect"
  ) +
  scale_x_continuous(
    limits = x_limits,
    sec.axis = sec_axis(
      trans = ~., 
      breaks = c(x_limits[1], x_limits[2]),
      labels = c("Depleted", "Enriched")
    )
  ) +
  theme_minimal() +
  theme(
    axis.text.x.top = element_text(face = "bold", color = "black"),
    axis.ticks.x.top = element_blank(),
    axis.title.x.top = element_blank()
  ) +
  scale_shape_manual(values = c(10,15,16,12,9,17,14,8,11))



