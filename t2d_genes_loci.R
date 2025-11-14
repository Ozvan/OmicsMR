setwd("/Users/ald73/Library/CloudStorage/OneDrive-UniversityofCambridge/ITG/T2DGGI_multiOmics")
library(data.table)
library(readxl)
library(dplyr)
library(GenomicRanges)
library(biomaRt)

# For causal genes, extract TSS gene in b37
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",
                host = "grch37.ensembl.org")

################################################################################
# ---------------------- Load tested and significant genes ---------------------
################################################################################
load("AllTestedGenes_eQTL_TissueSpecific.RData")
all.tested.dt <- rbindlist(lapply(names(res.alltested.alltissues), function(l1) {
                data.table(gene=res.alltested.alltissues[[l1]], source=gsub(" ", "_", l1))}), use.names = TRUE, fill=T)

load("AllTestedGenes_eQTL_pQTL_Blood_ENSG.RData")
all.tested.dt_list <- unlist(lapply(c("pQTL", "eQTL"), function(l1) {
  lapply(names(Tested.eQTL.pQTL.blood.ENSG[[l1]]), function(l2) {
    data.table(gene = Tested.eQTL.pQTL.blood.ENSG[[l1]][[l2]], source = paste(l1, l2, sep="_"))
  })}),recursive=F)

all.tested.dt <- rbind(all.tested.dt, rbindlist(all.tested.dt_list, use.names = TRUE, fill=T))
all.tested.dt <- all.tested.dt[, .(source = paste(unique(source), collapse = ",")), by = gene]

# Get HGNC symbol
all.tested.genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand"),
                  filters = "ensembl_gene_id",
                  values = unique(all.tested.dt$gene),
                  mart = mart)

all.tested.dt <- as.data.table(merge(all.tested.dt, all.tested.genes, by.x="gene", by.y="ensembl_gene_id"))

load("AllSigGenes_eQTL_TissueSpecific.RData")
dt <- rbindlist(lapply(names(res.sig.alltissues), function(l1) {
                data.table(gene=res.sig.alltissues[[l1]], source=gsub(" ", "_", l1))}), use.names = TRUE, fill=T)

load("AllSigGenes_eQTL_pQTL_rep_ENSG.RData")
dt_list <- unlist(lapply(c("pQTL", "eQTL"), function(l1) {
  lapply(names(all.significant.genes.eQTL.pQTL.ENSG[[l1]]), function(l2) {
    data.table(gene = all.significant.genes.eQTL.pQTL.ENSG[[l1]][[l2]], source = paste(l1, l2, sep="_"))
  })}),recursive=F)

dt <- rbind(dt, rbindlist(dt_list, use.names = TRUE, fill=T))
dt <- unique(dt[source!="eQTL_MetaAnalysis_FE"])
dt_unique <- dt[, .(source = paste(unique(source), collapse = ",")), by = gene]

# Get HGNC symbol
gene_pos <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"),
                  filters = "ensembl_gene_id",
                  values = unique(dt_unique$gene),
                  mart = mart)

gene_pos <- as.data.table(gene_pos)
gene_pos[, tss := ifelse(strand == 1, start_position, end_position)]
gene_pos <- gene_pos[chromosome_name %in% seq(1,22)]
gene_pos <- merge(gene_pos, dt_unique, by.x="ensembl_gene_id", by.y="gene")

# Define a 1Mb window (+/- 500kb) either side of TSS for each genes
window <- 500000
gr_causal <- GRanges(seqnames=gene_pos$chromosome_name,
                     ranges=IRanges(start=pmax(gene_pos$tss - window, 1), end=gene_pos$tss + window),
                     strand="*",
                     gene=gene_pos$ensembl_gene_id)

################################################################################
# -------------------- Load established T2D effector genes ---------------------
################################################################################
# For DM genes extract genes position in b37
dm.genes.all <- fread("dmgenes_unfiltered_GTExv10_21.08.25.txt")
dm.genes <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"),
                  filters = "external_gene_name",
                  values = unique(dm.genes.all$gene_id),
                  mart = mart)
dm.genes <- as.data.table(dm.genes)
dm.genes <- dm.genes[chromosome_name %in% seq(1,22)]
dm.genes <- merge(dm.genes, dm.genes.all, by.x="external_gene_name", by.y="gene_id")
dm.genes[, tested:=ifelse(external_gene_name %in% all.tested.dt$external_gene_name, 1, 0)]

dm_gr <- GRanges(seqnames=dm.genes$chromosome_name,
                 ranges=IRanges(start=dm.genes$start_position, end=dm.genes$end_position),
                 strand="*",
                 gene=dm.genes$ensembl_gene_id)

################################################################################
# -------------- Check overlap of ours and T2D effector genes ------------------
################################################################################
# For each causal gene look if any DM gene overlaps it
hits <- findOverlaps(gr_causal, dm_gr, ignore.strand = TRUE)
overlap.dt <- data.table(causal_gene=mcols(gr_causal)$gene[queryHits(hits)],
                         dm_gene=mcols(dm_gr)$gene[subjectHits(hits)])

overlap.dt <- unique(overlap.dt)
overlap.dt <- merge(overlap.dt, dt_unique, by.x = "causal_gene", by.y = "gene", all.x = TRUE)
overlap.dt
length(unique(overlap.dt$causal_gene))  # 582
overlap.dt[causal_gene==dm_gene]  # 134
matching.genes <- unique(overlap.dt[causal_gene==dm_gene, causal_gene])
overlap.dt[causal_gene %in% matching.genes]

overlap.dt[!(causal_gene %in% dm.genes$ensembl_gene_id)] 

# 1) New loci: putative with no overlaps within 1 Mb
put_idx_with_hit <- unique(queryHits(hits))
is_new <- !(seq_along(gr_causal) %in% put_idx_with_hit)
new_loci_df <- as_tibble(gene_pos[is_new, , drop = FALSE])

# 2) Overlapping map with equality / counts
ov_df <- tibble(
  put_index = queryHits(hits),
  cur_index = subjectHits(hits),
  put_gene  = mcols(gr_causal)$gene[queryHits(hits)],
  cur_gene  = mcols(dm_gr)$gene[subjectHits(hits)],
  put_chr   = as.character(seqnames(gr_causal))[queryHits(hits)],
  cur_chr   = as.character(seqnames(dm_gr))[subjectHits(hits)],
  put_start = start(gr_causal)[queryHits(hits)],
  put_end   = end(gr_causal)[queryHits(hits)]
) %>%
  mutate(equal_gene = put_gene == cur_gene)

################################################################################
# -------------- Create nice output tables ------------------
################################################################################
# Summary per putative gene
overlap_summary <- ov_df %>%
  group_by(put_gene, put_chr, put_start, put_end) %>%
  summarise(
    any_equal = any(equal_gene),
    n_curated_overlaps = n(),
    n_equal_overlaps = sum(equal_gene),
    n_non_equal_overlaps = n_curated_overlaps - n_equal_overlaps,
    non_equal_curated_genes = paste0(unique(cur_gene[!equal_gene]), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_curated_overlaps), put_chr, put_start)

# Global counts
total_putative <- length(gr_causal)
n_new <- nrow(new_loci_df)
n_overlap <- total_putative - n_new
n_equal_pairs <- if (nrow(ov_df) > 0) sum(ov_df$equal_gene) else 0
n_putative_with_equal <- if (nrow(overlap_summary) > 0) sum(overlap_summary$any_equal) else 0
n_putative_with_non_equal <- if (nrow(overlap_summary) > 0) sum(overlap_summary$n_non_equal_overlaps > 0) else 0

summary_counts <- tibble(total_putative = total_putative,
                         window_bp = window,
                         putative_with_any_overlap = n_overlap,
                         per_putative_with_any_overlap = n_overlap*100/total_putative,
                         putative_without_overlap_new_loci = n_new,
                         per_putative_without_overlap_new_loci = n_new*100/total_putative,
                         overlapping_pairs_equal_symbol = n_equal_pairs,
                         per_overlapping_pairs_equal_symbol = n_equal_pairs*100/total_putative,
                         putative_genes_with_at_least_one_equal_match = n_putative_with_equal,
                         per_putative_genes_with_at_least_one_equal_match = n_putative_with_equal*100/total_putative,
                         putative_genes_with_only_non_equal_matches = n_putative_with_non_equal,
                         per_putative_genes_with_only_non_equal_matches = n_putative_with_non_equal*100/total_putative)

################################################################################
# -------------- Add information of novel T2DGGI genetic loci ------------------
################################################################################
novel.t2d.signals <- fread("novel_t2dggi_snps.txt", header=FALSE)
snp_mart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp", GRCh = 37)
novel.t2d.signals <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
                     filters = "snp_filter",
                     values = novel.t2d.signals$V1,
                     mart = snp_mart)
novel.t2d.signals <- as.data.table(novel.t2d.signals)
novel.t2d.signals <- novel.t2d.signals[chr_name %in% seq(1,22)]

novel_snps_gr <- GRanges(seqnames = novel.t2d.signals$chr_name,
                   ranges = IRanges(start = novel.t2d.signals$chrom_start, end = novel.t2d.signals$chrom_start),
                   rsid = novel.t2d.signals$refsnp_id)

put_idx_with_hit <- unique(queryHits(hits))
is_new <- !(seq_along(gr_causal) %in% put_idx_with_hit)
new_loci_df <- as_tibble(gene_pos[is_new, , drop = FALSE])

# Check if new_loci_df overlap t2d.loci
new_loci_gr <- GRanges(
  seqnames = new_loci_df$chromosome_name,
  ranges = IRanges(
    start = new_loci_df$start_position,
    end = new_loci_df$end_position
  ),
  strand = ifelse(new_loci_df$strand == 1, "+", "-"),
  gene_id = new_loci_df$ensembl_gene_id,
  source = new_loci_df$source
)

# Find overlaps between new loci and SNPs
hits <- findOverlaps(novel_snps_gr, new_loci_gr)

# Genes that overlap at least one SNP
genes_with_overlap <- unique(new_loci_gr$gene_id[subjectHits(hits)])

# Genes in your loci that do not overlap any SNP
genes_without_overlap <- setdiff(unique(new_loci_df$ensembl_gene_id),
                                 genes_with_overlap)

# Number of unique non-overlapping genes
length(genes_without_overlap)

# Inspect overlapping pairs
overlaps <- data.frame(
  snp = mcols(snps_gr)$rsid[queryHits(hits)],
  chr = as.character(seqnames(snps_gr)[queryHits(hits)]),
  snp_pos = start(snps_gr)[queryHits(hits)],
  gene = new_loci_gr$gene_id[subjectHits(hits)],
  locus_chr = as.character(seqnames(new_loci_gr)[subjectHits(hits)]),
  locus_start = start(new_loci_gr)[subjectHits(hits)],
  locus_end = end(new_loci_gr)[subjectHits(hits)],
  source = new_loci_gr$source[subjectHits(hits)]
)

head(overlaps)
################################################################################
# -------------- Add tested and expressed information ------------------
################################################################################
# Add tissue of putative causal effect
overlap_summary_tissues <- merge(as.data.table(ov_df), dt_unique[, .(gene, source_put_gene=source)],  by.x="put_gene", by.y="gene")

# Add tested flag
overlap_summary_tissues <- merge(overlap_summary_tissues, dm.genes[, .(ensembl_gene_id, cur_gene_tested=tested)],  by.x="cur_gene", by.y="ensembl_gene_id")

# Add if curated gene expressed in tissue of putative causal effect
tissue_map <- list("Brain_Hypothalamus" = "Brain",
                  "Pancreas"           = "Pancreas",
                  "Pancreatic_Islets"  = "Islets",
                  "Liver"              = "Liver",
                  "Skeletal_muscle"    = "Muscle",
                  "Subcutaneous_Adipose" = "SAT",
                  "Visceral_Adipose" = "VAT",
                  "eQTL_EUR" = "WholeBlood",
                  "eQTL_MetaAnalysis" = "WholeBlood",
                  "pQTL_MetaAnalysis" = "WholeBlood",
                  "pQTL_EUR" = "WholeBlood",
                  "pQTL_EAS" = "WholeBlood",
                  "eQTL_AMR" = "WholeBlood",
                  "eQTL_AFR" = "WholeBlood")

# Make sure dm.genes is keyed on ensembl_gene_id for fast join
setkey(dm.genes, ensembl_gene_id)

# Add expressed flag
overlap_summary_tissues[, cur_gene_expressed := {
  if (is.null(tissue_map[[source_put_gene]])) {
    NA_character_
  } else {
    val <- dm.genes[cur_gene, get(tissue_map[[source_put_gene]])]
    ifelse(!is.na(val) && val > 0, "yes", "no")
  }
}, by = .(cur_gene, source_put_gene)]

equal.genes=overlap_summary_tissues[equal_gene==T, put_gene]
overlap_summary_tissues[equal_gene==F & !put_gene %in% equal.genes, uniqueN(put_gene)]
overlap_summary_tissues[equal_gene==F & !put_gene %in% equal.genes & cur_gene_tested==1, uniqueN(put_gene)]
overlap_summary_tissues[equal_gene==F & !put_gene %in% equal.genes & cur_gene_tested==1 & cur_gene_expressed=="yes", uniqueN(put_gene)]

################################################################################
# -------------- Look at distance to T2D signals ------------------
################################################################################
t2d.signals <- fread("t2dggi_risk_variants.txt", header=FALSE)

snp_mart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp", GRCh = 37)
t2d.signals <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
                     filters = "snp_filter",
                     values = t2d.signals$V1,
                     mart = snp_mart)
t2d.signals <- as.data.table(t2d.signals)
t2d.signals <- t2d.signals[chr_name %in% seq(1,22)]

snps_gr <- GRanges(seqnames = t2d.signals$chr_name,
                   ranges = IRanges(start = t2d.signals$chrom_start, end = t2d.signals$chrom_start),
                   rsid = t2d.signals$refsnp_id)

# --- Filter overlaps where causal != curated ---
non_equal_hits <- ov_df %>% filter(!equal_gene)
t2d_gr <- snps_gr

# --- Loop over each non-equal overlap to find closest T2D SNP in locus ---
closest_snp_list <- lapply(seq_len(nrow(non_equal_hits)), function(i){
  put_start <- non_equal_hits$put_start[i]
  put_end   <- non_equal_hits$put_end[i]
  put_chr   <- non_equal_hits$put_chr[i]
  
  # Subset T2D SNPs in the locus window
  locus_snps <- t2d_gr[seqnames(t2d_gr) == put_chr & start(t2d_gr) >= put_start & start(t2d_gr) <= put_end]
  
  if(length(locus_snps) == 0){
    return(data.table(
      put_gene = non_equal_hits$put_gene[i],
      cur_gene = non_equal_hits$cur_gene[i],
      closest_t2d_snp = NA,
      dist_to_put = NA,
      dist_to_cur = NA
    ))
  }
  
  # Find the closest SNP to the causal TSS
  distances_to_put <- abs(start(locus_snps) - non_equal_hits$put_start[i])
  closest_idx <- which.min(distances_to_put)
  closest_snp <- locus_snps[closest_idx]
  
  # Calculate distance to causal (put) and curated (cur) genes
  dist_to_put <- abs(start(closest_snp) - non_equal_hits$put_start[i])
  dist_to_cur <- abs(start(closest_snp) - start(dm_gr[dm_gr$gene == non_equal_hits$cur_gene[i]]))
  
  data.table(
    put_gene = non_equal_hits$put_gene[i],
    cur_gene = non_equal_hits$cur_gene[i],
    closest_t2d_snp = mcols(closest_snp)$rsid,
    dist_to_put = dist_to_put,
    dist_to_cur = dist_to_cur
  )
})

# Combine results
closest_snp_dt <- rbindlist(closest_snp_list)
head(closest_snp_dt)
closest_snp_dt[dist_to_put<dist_to_cur, .N]
closest_snp_dt[, uniqueN(put_gene)]
closest_snp_dt[dist_to_put<dist_to_cur, uniqueN(put_gene)]
closest_snp_dt[, uniqueN(cur_gene)]
closest_snp_dt[dist_to_put<dist_to_cur, uniqueN(cur_gene)]
