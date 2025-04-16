library(tidyverse)
library(data.table)
library(magrittr)
library(coloc)

# parameters #########
if(FALSE){
  wd <- #Path of working directory
  cohort_name <- #Basename of the cohort
  outcome_name <- 'T2DGGI' #Name of the outcome
  
  samplesize_qtl <- #Sample size of the QTL study
  casesize_t2d <- #Number of cases in the outcome GWAS
  samplesize_t2d<- #Total sample size in the outcome GWAS
  
  path_T2DGGI <- #Path to the T2DGGI summary statistics
  plink.path <- #Path to plink
  pwcoco.path <- #Path to pwcoco
  ref.bfile.path <- #Path to 1000 Genome bfiles
  cis_sumstats_dir <- #Path to cis- summary statistics
  tss.path <- #File with information on TSS
}

# NOTE: you have to modify the renaming code for the full QTL sumstats in the code below 
######################

######## parameters specified with commandArgs
# pQTL data
args <- commandArgs(trailingOnly=TRUE) 
cis_pqtl_path <- args[1] 
protname <- args[2] 
######## 

setwd(wd)

# get seqid from protname
seqid1 <- str_split(protname, pattern = '[_]')[[1]][1] # e.g., 3796
seqid2 <- str_split(protname, pattern = '[_]')[[1]][2] # e.g., 79
seqid <- paste0(seqid1, '_', seqid2)

shortname <- str_split(protname, pattern = '[_]')[[1]][3] # e.g., PCSK1
protname <- paste0(shortname, '.', seqid) # update the format (only needed for deCODE)
print(protname)

################
# To reduce computational time, proceed to pwcoco only if protein-disease association is nominally significant (pval < 0.05) in MR
################
# List all files in the directory 'output/or/'
or_files <- list.files('output/or/', full.names = TRUE) # Use full.names = TRUE to get the full file path

sum_df <- fread(paste0('output/or/', protname, '.or.txt'), fill = T)
sum_df %<>% rename(protein = id.exposure)

sum2 <- sum_df %>% filter(method == 'Inverse variance weighted' | method == 'Wald ratio') %>% group_by(protein) %>% arrange(pval) %>% filter(!duplicated(protein)) %>% ungroup()
sum2$shortname <- unlist(lapply(strsplit(sum2$protein, "[.]"), function(x)x[1]))  
sum2$seqid <- unlist(lapply(strsplit(sum2$protein, "[.]"), function(x)x[2]))  
sum2 %<>% dplyr::select(protein, seqid, shortname, method, nsnp, b, se, pval, lo_ci, up_ci) 

# proceed only if MR-IVW or Wald ratio's pval < 0.05 to reduce computational time
if(sum2$pval < 0.05){
  print(paste0('proceed to pwcoco for ', protname, '--', outcome_name))
  # create each protein's dir
  output_dir <- paste0(wd, '/output/pwcoco_res/', protname, '/')
  system(paste0('mkdir -p ', paste0(output_dir, protname, '--', outcome_name)))  # e.g., "/scratch/richards/satoshi.yoshiji/33.T2DGGI/pQTL_deCODE_plasma_EUR/output/pwcoco_res/ANGPTL4.3796_79/ANGPTL4.3796_79--T2DGGI"
  
  # extract the IV with the lowerst pval (lead cis-pQTL) from harmonized data
  harmonize_df <-  fread(paste0('output/harmonized/', protname, '.harmonized.txt'), fill = T)
  lead_cisqtl_rsid <- harmonize_df %>% arrange(pval.exposure) %>% filter(row_number() == 1)  %>% pull(SNP) # IV with the lowest pval

  ########### (ii) corresponding cis_sumstats ###############
    cis_sumstats_path <- system(paste0('ls ', cis_sumstats_dir, '/*' , seqid, '*.gz'), intern = T) 

    #check all
    print(protname)
    print(cis_sumstats_path)
    
    
    # NEED MODIFICATION ACCORDING TO EACH SUMSTATS' FORMAT
    cis_sumstats <- fread(cis_pqtl_path, fill = T) 
    cis_sumstats %<>% transmute(chromosome = Chrom, 
                       base_pair_location =Pos,
                       p_value = Pval, 
                       beta = Beta,
                       standard_error = SE, 
                       effect_allele_frequency = ImpMAF, 
                       # n, 
                       rs_id = rsids, 
                       effect_allele = effectAllele, 
                       other_allele = otherAllele)
    
    cis_sumstats %<>% mutate(n = samplesize_qtl)
    cis_sumstats$effect_allele <- str_to_upper(cis_sumstats$effect_allele) # for fenland
    cis_sumstats$other_allele <- str_to_upper(cis_sumstats$other_allele) # for fenland
    
    ##############
    # define the region for colocalization based on TSS plus minus 1 mb
    ##############
    tss_df <- readxl::read_excel(tss.path, sheet = 'ST6', skip = 1)
    tss_chr <- tss_df %>% filter(SeqId == seqid) %>% pull('Chromosome') %>% as.numeric()
    tss_pos <- tss_df %>% filter(SeqId == seqid) %>% pull('TSS')  %>% as.numeric()  
    # lead_cisqtl_pos <- cis_sumstats %>% filter(rs_id ==  lead_cisqtl_rsid) %>% pull(base_pair_location)
    cis_sumstats_sel <- cis_sumstats %>%
      filter(chromosome == paste0('chr', tss_chr)) %>%
      filter(base_pair_location >= tss_pos - 500000) %>%
      filter(base_pair_location <= tss_pos + 500000)
    
    ##############
    # in case where you define cis_sumstats based on the lead variant, use this
    ##############
    # # now, extract the lead cis-pQTL position
    # lead_cisqtl_chr <- cis_sumstats %>% filter(rs_id ==  lead_cisqtl_rsid) %>% pull(chromosome)
    # # lead_cisqtl_pos <- cis_sumstats %>% filter(rs_id ==  lead_cisqtl_rsid) %>% pull(base_pair_location)
    # cis_sumstats_sel <- cis_sumstats %>% 
    #   filter(chromosome == lead_cisqtl_chr) %>%
    #   filter(base_pair_location >= lead_cisqtl_pos-500000) %>%
    #   filter(base_pair_location <= lead_cisqtl_pos+500000)
    
    # give varbeta
    cis_sumstats_sel <- cis_sumstats_sel %>% mutate(varbeta = standard_error^2)
    
    # remove variants w/o rsid
    cis_sumstats_rsid <- cis_sumstats_sel %>%
      filter(!is.na(rs_id))
    # remove duplicated snp
    cis_sumstats_nodup <- cis_sumstats_rsid %>%
      dplyr::arrange(p_value) %>%
      dplyr::filter(!duplicated(rs_id))
    
    # rename 
    cis_sumstats_nodup %<>%
      transmute(
        chr = chromosome,
        position = base_pair_location,
        beta = beta,
        varbeta = varbeta,
        SE = standard_error,
        ALT = effect_allele,
        REF = other_allele,
        MAF = effect_allele_frequency,
        N = n,
        pvalues = p_value,
        snp = rs_id
      ) 
    
    # rename to cis_sumstats_nodup and add id column
    cis_sumstats_nodup$id <- paste('chr', cis_sumstats_nodup$chr, cis_sumstats_nodup$position, cis_sumstats_nodup$REF, cis_sumstats_nodup$ALT, sep = ":")
    
    # drop_na
    #cis_sumstats_nodup <- drop_na(cis_sumstats_nodup)
    cis_sumstats_nodup <- cis_sumstats_nodup[complete.cases(cis_sumstats_nodup), ]
    
    # deal w/ extreme MAF and pval
    cis_sumstats_nodup$MAF[cis_sumstats_nodup$MAF == 0] <- 1e-4
    cis_sumstats_nodup$MAF[cis_sumstats_nodup$MAF == 1] <- 0.9999
    
    cis_sumstats_nodup$pvalues[cis_sumstats_nodup$pvalues == 0] <- 1e-300
    cis_sumstats_nodup$pvalues[cis_sumstats_nodup$pvalues == 1] <- 0.9999
    
    # # make it a list
    # npnt_list <- as.list(npnt_nodup)
    # npnt_list$type <- 'quant'
    # str(npnt_list)
    #
    # # check dataset
    # check_dataset(npnt_list)
    
    #####################
    # outcome GWAS
    #####################
    outcome <- fread(path_T2DGGI, fill = T)
    outcome_rename0 <- outcome %>%
      dplyr::rename(
        chr = CHR_HG38,
        position = POS_HG38,
        beta = Effect,
        SE = StdErr,
        ALT = Allele2,
        REF = Allele1,
        pvalues = 'P-value',
        MAF = Freq1,
        #  samplesize = SS,
        snp = RSID)
    
    # keep SNPs that are present in cis_sumstats
    outcome_rename <- outcome_rename0 %>% filter(snp %in% unlist(cis_sumstats_nodup$snp))  # to reduce computational time, restrict SNPs here
    
    outcome_rename %<>% mutate(id = paste0('chr', chr, ':', position, ':', REF, ':', ALT)) #ALT = effect_allele, REF = other_allele,
    
    outcome_sel <- outcome_rename %>% dplyr::select(id, chr, position, beta, SE, ALT, REF, pvalues, MAF, snp)
    
    # give varbeta
    outcome_sel$varbeta <- (outcome_sel$SE)^2
    #outcome_sel %<>% dplyr::select(-SE)
    
    # remove variants w/o rsid
    outcome_rsid <- outcome_sel %>%
      filter(!is.na(snp))
    
    # remove duplicated id
    outcome_nodup <- outcome_rsid %>%
      dplyr::arrange(pvalues) %>%
      dplyr::filter(!duplicated(snp))
    
    # remove NA
    #outcome_nodup <- drop_na(outcome_nodup)
    
    # deal w/ extreme MAF and pval to avoid errors
    outcome_nodup$MAF[outcome_nodup$MAF == 0] <- 1e-4
    outcome_nodup$MAF[outcome_nodup$MAF == 1] <- 0.9999
    
    outcome_nodup$pvalues[outcome_nodup$pvalues == 0] <- 1e-300
    outcome_nodup$pvalues[outcome_nodup$pvalues == 1] <- 0.9999
    
    ####################
    # keep only common SNPs
    # inner_join will not only keep rows with common SNPs but make their order same in both datasets
    ####################
    print('working on common')
    
    common <- inner_join(cis_sumstats_nodup, outcome_nodup, by = 'snp', suffix = c('.cis_sumstats', '.outcome'))
    common %<>% rename(ALT = ALT.cis_sumstats, REF = REF.cis_sumstats) %>% select(-ALT.outcome) %>% select(-REF.outcome)
    
    
    # If there is no common SNPs, write error message out and skip the process. otherwise the entire process will stop.
    if(nrow(common) == 0){
      #resname <- paste0(protname, '.', outcome_name, '.coloc.tsv')
      #write.table('no common snp', file = resname, quote = F, row.names = F, sep ='\t')
      print('no common snp')
    } else {
      print('common snps exist')
      # pQTL
      cis_sumstats_common <- common %>% dplyr::select(chr.cis_sumstats, position.cis_sumstats, ALT, REF, beta.cis_sumstats,  pvalues.cis_sumstats, MAF.cis_sumstats, varbeta.cis_sumstats, SE.cis_sumstats, snp) %>%
        dplyr::rename(chr = chr.cis_sumstats, position = position.cis_sumstats, ALT = ALT, REF = REF, beta = beta.cis_sumstats, pvalues = pvalues.cis_sumstats, MAF = MAF.cis_sumstats, varbeta = varbeta.cis_sumstats, SE = SE.cis_sumstats)
      
      # make sure to obtain chr and pos derives from cis_sumstats (not outcome, which can be based on GRCh37/38)
      outcome_common <- common %>% dplyr::select(chr.cis_sumstats, position.cis_sumstats, ALT, REF, beta.outcome, pvalues.outcome, MAF.outcome, varbeta.outcome, SE.outcome, snp) %>%
        dplyr::rename(chr = chr.cis_sumstats, position = position.cis_sumstats, ALT = ALT, REF = REF, beta = beta.outcome, pvalues = pvalues.outcome, MAF = MAF.outcome, varbeta = varbeta.outcome, SE = SE.outcome)
      
      # run pwcoco only if there are SNPs with pval <0.05 in both traits to reduce computational time
      if(min(cis_sumstats_common$pvalues) < 0.05 & min(outcome_common$pvalues) < 0.05){
        ######
        # run pwcoco
        ######
        sumstats1 <- cis_sumstats_common
        sumstats2 <- outcome_common
        
        ##################
        # pQTL
        ##################
        # sumstats1_list <- as.list(sumstats1)
        # sumstats1_list$type <- 'quant'
        # sumstats1_list$N <- cispqtl_lowestpval$n
        # sumstats1_list$MAF <- sumstats1$MAF
        # #sumstats1_list$sdY <- 1
        # str(sumstats1_list)
        # # check dataset
        # check_dataset(sumstats1_list)
        # plot_dataset(sumstats1_list)
        # 
        #pwcoco
        sumstats1_pwcoco <- sumstats1 %>% select(snp, ALT, REF, MAF, beta, SE, pvalues)
        sumstats1_pwcoco %<>% mutate(n = samplesize_qtl)
        sumstats1_name <- paste0(output_dir, protname, '--', outcome_name, '/', protname, '--', outcome_name, '--sumstats1.tsv')
        system(paste0('mkdir -p ', paste0(output_dir, protname, '--', outcome_name)))
        write_tsv(sumstats1_pwcoco, file = sumstats1_name)

        #pwcoco
        sumstats2_pwcoco <- sumstats2 %>% select(snp, ALT, REF, MAF, beta, SE, pvalues)
        sumstats2_pwcoco %<>% mutate(n = as.integer(samplesize_t2d)) %>% mutate(case = as.integer(casesize_t2d))
        sumstats2_name <-paste0(output_dir, protname, '--', outcome_name, '/', protname, '--', outcome_name, '--sumstats2.tsv')
        write_tsv(sumstats2_pwcoco, file = sumstats2_name)
      }
      

      # Now both dataset 1 and 2 are ready, so let's perform pwcoco (instead of coloc)
      # pwcoco
      system(paste0(pwcoco.path, ' --bfile ', ref.bfile.path  , ' --sum_stats1 ',  sumstats1_name,  ' --sum_stats2 ', sumstats2_name, ' --top_snp 5 --maf 0.01 --out ', output_dir, protname, '--', outcome_name, '/', protname, '--', outcome_name))
      print('pwcoco completed')
      
      
      # if you want to run coloc instead of pwcoco, use below 
      # pwcoco_res <- coloc.abf(
      #   dataset1 = sumstats1_list,
      #   dataset2 = sumstats2_list
      # )
      
      # # write results
      # system(paste0('mkdir -p ', protname))
      # resname <- paste0(protname, '/', protname, '-', outcome_name, '-coloc.tsv')
      # write_tsv(as.data.frame(pwcoco_res$summary) %>% rownames_to_column(), file = resname)
    } # close    if(nrow(common) == 0){

    
} else { print(paste0('skip pwcoco for ', protname, '--', outcome_name))}    # closes if(sum2$pval < 0.05){
