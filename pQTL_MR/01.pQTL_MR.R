library(TwoSampleMR)
library(magrittr)
library(tidyverse)
library(data.table) 
library(MRPRESSO)

######## parameters to fill:
if(FALSE){
  wd <- #Path to data
  cohort_name <- #Basename of the cohort
  path_T2DGGI <- #File and path to the T2DGGI summary statistics
  samplesize_qtl <- #Sample size of the QTL dataset
  samplesize_t2d <- #Sample size of the outcome (T2D) GWAS meta-analysis
  plink.path <- #Path to plink
  ref.bfile.path <- #Path to 1000 Genomes project
  pval.thr <- #Significance threshold
  r2_thr <- #r2 value for pruning
  kb_thr <- #window size for pruning
}
###################

######## parameters specified with commandArgs
args <- commandArgs(trailingOnly=TRUE) 
qtl.path <- args[1] # this refers to conditionally independent and significant SNPs defined by the original consortium. 
qtl.name <- args[2] 
###################

######## another parameter specified by seqid, which can be extracted from qtl.name (e.g., seqid = 620_13). This has to be modified per cohort
seqid <- stringr::str_extract(qtl.name, "\\d+_\\d+")
cis_sumstats_path <- system(paste0('ls /scratch/richards/satoshi.yoshiji/11.pQTL/10.deCODE_cis_full_sumstats/sumstats/*', seqid, '*.cis.txt.gz'), intern = T)

system(paste0('mkdir -p ', wd))
setwd(wd)
output_dir <- paste0(wd, "/output/")
system(paste0('mkdir -p ', output_dir, 'harmonized/'))
system(paste0('mkdir -p ', output_dir, 'mr/'))
system(paste0('mkdir -p ', output_dir, 'or/'))
system(paste0('mkdir -p ', output_dir, 'hetero/'))
system(paste0('mkdir -p ', output_dir, 'pleio/'))
system(paste0('mkdir -p ', output_dir, 'steiger/'))
system(paste0('mkdir -p ', output_dir, 'mrpresso/'))


# read the outcome GWAS
t2d <- fread(path_T2DGGI, fill = T)

# exposure
exp_dat <- read_exposure_data(filename = qtl.path,
                              sep='\t',
                              snp_col='variant',
                              beta_col = 'beta_unadj',
                              se_col = 'se',
                              effect_allele_col = 'Amin',
                              other_allele_col = 'Amaj',
                              eaf_col = 'MAF',
                              pval_col = 'pval'#,
                              #samplesize_col = "N"
                              )
exp_dat %<>% mutate(rsid = SNP) %>% mutate(samplesize.exposure = samplesize_qtl)

##################
# Below, we try to find a proxy 
##################

# clump SNPs
system(paste0("mkdir -p ", wd, "/Clumping/", qtl.name)) # e.g., /scratch/richards/satoshi.yoshiji/33.T2DGGI/pQTL_deCODE_plasma_EUR/Clumping/PAM.5620_13
write_tsv(exp_dat, file = paste0(wd, "/Clumping/", qtl.name, '/', qtl.name, '.tsv'))

output.qtl <- data.table()
# The firt arguments of the functions are for the clumping, and the last ones for the data
# Clump the variants

system(paste0(plink.path, " --bfile ", ref.bfile.path, " --clump-p1 ", pval.thr, " --clump-r2 ", r2_thr, " --clump-kb ", kb_thr, " --clump ",  paste0(wd, "/Clumping/", qtl.name, '/', qtl.name, '.tsv'), " --clump-snp-field rsid --clump-field pval.exposure -out ", wd, "/Clumping/", qtl.name, '/', qtl.name))
  #Import the clumping results
  #Need to check if there are variants left after clumping, not the case for example if only ine significant variant but not in reference data
  if(!file.exists(paste0(wd, "/Clumping/", qtl.name, '/', qtl.name, '.clumped'))) return(NA) # e.g., /scratch/richards/satoshi.yoshiji/33.T2DGGI/pQTL_deCODE_plasma_EUR/Clumping/PAM.5620_13/PAM.5620_13.clumped

clumping.res <- fread(paste0(wd, "/Clumping/", qtl.name, '/', qtl.name, ".clumped"))
clumped.ivs <- clumping.res$SNP

# cis sumstats
qtl <- fread(cis_sumstats_path, fill = T) 
qtl %<>% transmute(chr=Chrom, 
                   pos=Pos,
                   pval=Pval, 
                   beta=Beta,
                   se=SE, 
                   eaf= ImpMAF,
                   rsid= rsids, 
                   ea=effectAllele, 
                   nea=otherAllele)

#Import the clumping results
#Need to check if there are variants left after clumping, not the case for example if only ine significant variant but not in reference data
if(nrow(qtl) == 0) return(NA)

#Compute F-stat
qtl$snp_fstat <- (qtl$beta)^2/(qtl$se)^2

#Only keep the variants that are clumped
qtl <- qtl[, iv:=ifelse(rsid %in% clumped.ivs & snp_fstat>=10, TRUE, FALSE)]

##################
# proxy search function
####################
#Check which clumped IVs are not available in T2D of matching ancestry
need.proxy <- subset(qtl, iv)$rsid[which(!(subset(qtl, iv)$rsid %in% t2d$RSID))]

## Find proxies if IV not included in T2D sumstats
if(length(need.proxy)>0){
  proxies <- NULL
  iv=need.proxy[1]
  for (iv in need.proxy){
    cluster.snps <- unlist(strsplit(gsub("\\(1\\)", "", clumping.res[SNP==iv, SP2]), split=","))
    
    # Check LD between iv and potential.proxy
    ld.matrix <- ieugwasr::ld_matrix(c(iv, cluster.snps),
                                     plink_bin = genetics.binaRies::get_plink_binary(),
                                     bfile = ref.bfile.path)
    #Take the corresponding iv and remove the correlation with itself
    ld.vector <- ld.matrix[grep(rownames(ld.matrix), pattern = iv),-grep(rownames(ld.matrix), pattern = iv)]
    rm(ld.matrix)
    
    for (potential.proxy in sub("_.*", "", names(which((ld.vector**2)>0.8)))){
      # Check if potential.proxy is in T2D data
      if(nrow(t2d[rs_id==potential.proxy])!=0) break
    }
    proxies <- c(proxies, potential.proxy)
  }
  
  ## Update IVs in sig.qtl
  qtl[rsid %in% need.proxy, iv:=FALSE]
  qtl[rsid %in% proxies, iv:=TRUE]
  
}

output.qtl <- rbind(output.qtl, qtl[iv==TRUE]) # this will be the proxy to be used
exposure_dat <-format_data(output.qtl, 
                           type="exposure", snp_col="rsid", beta_col="beta", se_col = "se", eaf_col = "eaf", 
                           effect_allele_col = "ea", other_allele_col = "nea", pval_col = "pval", chr_col = "chr", pos_col = "pos") # there is no sample size column
exposure_dat$exposure <- qtl.name
exposure_dat$id.exposure <- qtl.name
exposure_dat$samplesize.exposure <- samplesize_qtl

#######################
# common 
# downstream procedures are the same from MR with or without the use of proxies
#######################

# format the outcome
t2d$Allele1 <- str_to_upper(t2d$Allele1) 
t2d$Allele2 <-  str_to_upper(t2d$Allele2)

formatted_outcome <- format_data(t2d, 
                                 snps = output.qtl$rsid, # the function will extract SNPs only in output.qtl$rsid. This will speed up the process
                                 type="outcome", snp_col="RSID", beta_col="Effect", se_col = "StdErr", eaf_col = "Freq1", 
                                 effect_allele_col = "Allele1", other_allele_col = "Allele2", pval_col = "P-value", chr_col = "CHR_HG38", pos_col = "POS_HG38") # there is no sample size column
formatted_outcome$id.outcome <- 'outcome'
formatted_outcome$samplesize.outcome <- samplesize_t2d 

# harmonize
exp_dat_outcome <-harmonise_data(exposure_dat=exposure_dat, outcome_dat=formatted_outcome)
exp_data_outcome_name <- paste0(output_dir,"harmonized/", qtl.name, ".harmonized.txt")
write.table(exp_dat_outcome, file=exp_data_outcome_name, sep = '\t', quote = F, row.names = F)

# mr result
mr_results <- mr(exp_dat_outcome)
mr_name <- paste0(output_dir, "mr/", qtl.name, ".mr.txt")
write.table(mr_results, file=mr_name, sep = '\t', quote = F)

# odds ratio
OR <- generate_odds_ratios(mr_results)
OR_name <- paste0(output_dir, "or/", qtl.name, ".or.txt")
write.table(OR, file=OR_name, sep = '\t', quote = F, row.names = F)

# # scatter plot
# pdf_name <- paste0(output_dir, "pdf/", qtl.name, ".scatter.pdf")
# pdf(pdf_name)
# mr_scatter_plot(mr_results, exp_dat_outcome)[[1]] + xlim(0,1)
# dev.off()

# horizontal pleiotropy
tryCatch({ 
  pleio_res <- mr_pleiotropy_test(exp_dat_outcome)
  pleio_name <- paste0(output_dir,"pleio/", qtl.name, ".pleio.txt")
  write.table(pleio_res, file=pleio_name, sep = '\t', quote = F, row.names = F)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

# hetero test
tryCatch({ 
  hetero_res <- mr_heterogeneity(exp_dat_outcome)
  hetero_res$isquared <- abs(100*(hetero_res$Q - hetero_res$Q_df)/hetero_res$Q)  # I2 = 100%×(Q - df)/
  hetero_name <- paste0(output_dir,"hetero/", qtl.name, ".hetero.txt")
  write.table(hetero_res, file=hetero_name, sep = '\t', quote = F, row.names = F)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

# steiger
tryCatch({ 
  steiger <- directionality_test(exp_dat_outcome)
  steiger_name <- paste0(output_dir, "steiger/", qtl.name, ".steiger.txt")
  write.table(steiger, file=steiger_name, sep = '\t', quote = F, row.names = F)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

# MR-PRESSO
tryCatch({ 
  mrpresso_res <- mr_presso(BetaExposure = "beta.exposure",
                            BetaOutcome = "beta.outcome",
                            SdOutcome = "se.outcome",
                            SdExposure = "se.exposure",
                            OUTLIERtest = TRUE,
                            DISTORTIONtest = TRUE,
                            data = exp_dat_outcome,
                            NbDistribution = 1000,
                            SignifThreshold = 0.05)
  mrpresso_name <- paste0(output_dir, "mrpresso/", qtl.name, ".mrpresso.txt")
  
  # add global_rss, global_pval, distortion_indices, distortion_coef, distortion_pval
  mrpresso_df <- as.data.frame(mrpresso_res$`Main MR results`)
  mrpresso_df %<>% mutate(#global_rss = mrpresso_res$`MR-PRESSO results`$`Global Test`$RSSobs,
    global_pval =  mrpresso_res$`MR-PRESSO results`$`Global Test`$Pvalue,
    #distorition_indices = mrpresso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`,
    #distortion_coef = mrpresso_res$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`,
    distortion_pval = mrpresso_res$`MR-PRESSO results`$`Distortion Test`$Pvalue)
  
  mrpresso_df %<>% mutate(
    mrpresso_distortion_indices = ifelse(is.null(mrpresso_distortion_indices), NA, mrpresso_distortion_indices),
    mrpresso_distortion_coef = ifelse(is.null(mrpresso_distortion_coef), NA, mrpresso_distortion_coef),
    mrpresso_distortion_pval = ifelse(is.null(mrpresso_distortion_pval), NA, mrpresso_distortion_coef)  # Note: this line replaces NA with mrpresso_distortion_coef, which might be a mistake. If you meant to keep NA, use mrpresso_distortion_pval in the replacement.
  )
  
  write_tsv(mrpresso_df, file =  mrpresso_name)      
  
  
  #   # additional information
  #   mrpresso_global_rss <- mrpresso_res$`MR-PRESSO results`$`Global Test`$RSSobs
  #   mrpresso_global_pval <- mrpresso_res$`MR-PRESSO results`$`Global Test`$Pvalue
  #   mrpresso_distortion_indices <-  mrpresso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
  #   mrpresso_distortion_coef <-  mrpresso_res$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`
  #   mrpresso_distortion_pval <-  mrpresso_res$`MR-PRESSO results`$`Distortion Test`$Pvalue
  # 
  #   mrpresso_add_res <- data.frame(mrpresso_global_rss = mrpresso_global_rss,
  #                                  mrpresso_global_pval = mrpresso_global_pval)
  # mrpresso_add_res %<>% mutate(
  #   mrpresso_distortion_indices = ifelse(is.null(mrpresso_distortion_indices), NA, mrpresso_distortion_indices),
  #   mrpresso_distortion_coef = ifelse(is.null(mrpresso_distortion_coef), NA, mrpresso_distortion_coef),
  #   mrpresso_distortion_pval = ifelse(is.null(mrpresso_distortion_pval), NA, mrpresso_distortion_coef)  # Note: this line replaces NA with mrpresso_distortion_coef, which might be a mistake. If you meant to keep NA, use mrpresso_distortion_pval in the replacement.
  # )
  # 
  #   mrpresso_add_name <- paste0(output_dir, "mrpresso_add/", qtl.name, ".mrpresso_add.txt")
  #   write.table(mrpresso_add_res, file=mrpresso_add_name, sep = '\t', quote = F, row.names = F)
  
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

