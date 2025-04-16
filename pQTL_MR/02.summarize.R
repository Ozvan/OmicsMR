library(data.table)
library(tidyverse)
library(magrittr)
library(ggrepel)

# parameters #########
if(FALSE){
  wd <- #Path to working directory
  cohort_name <- #Basename of the cohort considered
}
######################
setwd(wd)
# List all files in the directory 'output/or/'
or_files <- list.files('output/or/', full.names = TRUE) # Use full.names = TRUE to get the full file path

# Initialize an empty data table for accumulating results
sum_df <- data.table()
# Loop through each file and combine their content into sum_df
for(ith_file in or_files) {
  # Read the current file
  ith_res <- fread(ith_file, fill = T)
  # Combine the data from ith_res into sum_df
  sum_df <- rbindlist(list(sum_df, ith_res), use.names = TRUE, fill = TRUE)
}
sum_df %<>% rename(protein = id.exposure)

#################
# get causal estimates
#################
# IVW or wald
sum2 <- sum_df %>% filter(method == 'Inverse variance weighted' | method == 'Wald ratio') %>% group_by(protein) %>% arrange(pval) %>% filter(!duplicated(protein)) %>% ungroup()
sum2$shortname <- unlist(lapply(strsplit(sum2$protein, "[.]"), function(x)x[1]))  
sum2$seqid <- unlist(lapply(strsplit(sum2$protein, "[.]"), function(x)x[2]))  
sum2 %<>% dplyr::select(protein, seqid, shortname, method, nsnp, b, se, pval, lo_ci, up_ci) 

# weighted median
sum2_median <- sum_df %>% filter(method == 'Weighted median') %>% group_by(protein) %>% arrange(pval) %>% filter(!duplicated(protein)) %>% ungroup()
sum2_median %<>% dplyr::select(protein, b, se, pval, lo_ci, up_ci) %>% dplyr::rename(b.median = b, se.median = se, lo_ci.median = lo_ci, up_ci.median = up_ci, pval.median = pval)

# egger
sum2_egger <- sum_df %>% filter(method == 'MR Egger') %>% group_by(protein) %>% arrange(pval) %>% filter(!duplicated(protein)) %>% ungroup()
sum2_egger %<>% dplyr::select(protein, b, se, pval, lo_ci, up_ci) %>% dplyr::rename(b.egger = b, se.egger = se, lo_ci.egger = lo_ci, up_ci.egger = up_ci, pval.egger = pval)

#################
# get MR-PRESSO results
#################
mrpresso_files <- list.files('output/mrpresso/', full.names = TRUE) # Use full.names = TRUE to get the full file path

# Initialize an empty data table for accumulating results
mrpresso_df <- data.table()

# Loop through each file
for(ith_file in mrpresso_files) {
  # Check if the file is not empty by ensuring its size is greater than 0
  if(file.info(ith_file)$size > 0) {
    # Read the current file if it is not empty
    ith_res <- vroom::vroom(ith_file, delim = '\t')
    
    # exttract protein name (with seqid) from the file name e.g., SMOC1.13118_5
    ith_protname <- sub("\\.mrpresso\\.txt$", "", basename(ith_file))
    
    # reorganize
    df <- cbind(data.frame(protein = ith_protname),
                ith_res[1,] %>% transmute(mrpresso.causal_estimate.raw = `Causal Estimate`, mrpresso.sd.raw = Sd, mrpresso.tstat.raw = `T-stat`, mrpresso.pval.raw = `P-value`, mrpresso.global_pval = global_pval),
                ith_res[2,] %>% transmute(mrpresso.causal_estimate.corrected = `Causal Estimate`, mrpresso.sd.corrected = Sd, mrpresso.tstat.corrected = `T-stat`, mrpresso.pval.corrected = `P-value`))
    ith_res %<>% mutate(protein = ith_protname, .before = `MR Analysis`) %>% transmute(protein, method.mrpresso  = `MR Analysis`, beta.mrpresso = `Causal Estimate`, sd.mrpresso = Sd, tstat.mrpresso = `T-stat`, pval.mrpresso = `P-value`, global_pval.mrpresso = global_pval)
    
     }
  # rbind
  mrpresso_df <- rbind(mrpresso_df, df)
}
mrpresso2 <- mrpresso_df

#################
# get hetero results
#################
hetero_files <- list.files('output/hetero/', full.names = TRUE) # Use full.names = TRUE to get the full file path

# Initialize an empty data table for accumulating results
hetero_df <- data.table()

# Loop through each file
for(ith_file in hetero_files) {
  # Check if the file is not empty by ensuring its size is greater than 0
  if(file.info(ith_file)$size > 0) {
    # Read the current file if it is not empty
    ith_res <- vroom::vroom(ith_file, delim = '\t')
    
    # Combine the data from ith_res into hetero_df
    hetero_df <- rbindlist(list(hetero_df, ith_res), use.names = TRUE, fill = TRUE)
  }
  # If the file is empty, the loop simply moves on to the next file
}
hetero2 <- hetero_df %>% filter(!is.na(Q_pval))
hetero2 %<>% rename(protein = id.exposure)
hetero2 %<>% dplyr::select(protein, method, Q, Q_df, Q_pval, isquared) 

#################
# get pleio results
#################
pleio_files <- list.files('output/pleio/', full.names = TRUE) # Use full.names = TRUE to get the full file path

# Initialize an empty data table for accumulating results
pleio_df <- data.table()

# Loop through each file
for(ith_file in pleio_files) {
  # Check if the file is not empty by ensuring its size is greater than 0
  if(file.info(ith_file)$size > 0) {
    # Read the current file if it is not empty
    ith_res <- vroom::vroom(ith_file, delim = '\t')
    
    # Combine the data from ith_res into pleio_df
    pleio_df <- rbindlist(list(pleio_df, ith_res), use.names = TRUE, fill = TRUE)
  }
  # If the file is empty, the loop simply moves on to the next file
}

pleio2 <- pleio_df %>% filter(!is.na(pval))
pleio2 %<>% rename(protein = id.exposure)
pleio2 %<>% dplyr::select(protein, egger_intercept, se, pval) %>% dplyr::rename(se.egger_intercept = se, pval.egger_intercept = pval)

#################
# get Steiger results
#################
steiger_files <- list.files('output/steiger/', full.names = TRUE) # Use full.names = TRUE to get the full file path

# Initialize an empty data table for accumulating results
steiger_df <- data.table()

# Loop through each file
for(ith_file in steiger_files) {
  # Check if the file is not empty by ensuring its size is greater than 0
  if(file.info(ith_file)$size > 0) {
    # Read the current file if it is not empty
    ith_res <- vroom::vroom(ith_file, delim = '\t')
    
    # Combine the data from ith_res into steiger_df
    steiger_df <- rbindlist(list(steiger_df, ith_res), use.names = TRUE, fill = TRUE)
  }
  # If the file is empty, the loop simply moves on to the next file
}

steiger2 <- steiger_df %>% filter(!is.na(steiger_pval))
steiger2 %<>% rename(protein = id.exposure)
steiger2 %<>% dplyr::select(protein, correct_causal_direction, steiger_pval) 

##############
# now start joining data frames
##############
# joining all IVW data first
if(is_empty(mrpresso2)){
  join <- left_join(sum2, hetero2 %>% filter(method == 'Inverse variance weighted') %>% select(-method), by='protein') %>%left_join(., pleio2, by = 'protein') # nested left join
  join %<>% mutate(mrpresso.causal_estimate.raw = NA, mrpresso.causal_estimate.corrected = NA)
} else {
  join <- left_join(sum2, hetero2 %>% filter(method == 'Inverse variance weighted') %>% select(-method), by='protein') %>% left_join(., mrpresso2, by = 'protein') %>%left_join(., pleio2, by = 'protein') # nested left join
  
}

# then join median,  egger and steiger results
join2 <- left_join(join, sum2_median, by='protein') %>% left_join(., sum2_egger, by = 'protein') %>% left_join(., steiger2, by = 'protein')  # nested left join

# remove duplicated (to ensure there is no duplicated protein)
join3 <- join2 %>% group_by(protein) %>% arrange(pval) %>% filter(!duplicated(protein))

########################
# annotate pass or not
########################

# Significance
join4 <- join3 %>% mutate(adjusted_pval = p.adjust(pval))

join4 %<>% 
  mutate(Significance = case_when(
    adjusted_pval < 0.05 ~ 'FDR',
    pval < 0.05 ~ 'Nominal',
    TRUE ~ 'No'
  ))

join4$Significance <- factor(join4$Significance, levels = c('No', 'Nominal', 'FDR'))

# count
table(join4$Significance)
hist(join4$pval)

# heterogeneity
# note that pass is defined as the followings:
# (i) Q_pval > 0.05
# (ii) Q_pval < 0.05 but weighted median, weighted mode, and MR-egger showed directionally consistent and FDR-significant results (0.05/4907)

join4 %<>% mutate(Heterogeneity_test = case_when(
                                                 isquared >= 50 ~ 'Fail',
                                                 isquared < 50 ~ 'Pass',
                                                 is.na(isquared) ~ 'NA'))
# count 
table(join4$Heterogeneity_test)

# Outlier robust estimate
join4 %<>% mutate(Outlier_robust_estimate_test = case_when(
  b > 0 & 
    (!is.na(b.median) & b.median > 0 | is.na(b.median)) & 
    (!is.na(b.egger) & b.egger > 0 | is.na(b.egger)) & 
    (!is.na(mrpresso.causal_estimate.raw) & mrpresso.causal_estimate.raw > 0 | is.na(mrpresso.causal_estimate.raw)) & 
    (!is.na(mrpresso.causal_estimate.corrected) & mrpresso.causal_estimate.corrected > 0 | is.na(mrpresso.causal_estimate.corrected)) ~ 'Pass',
  
  b < 0 & 
    (!is.na(b.median) & b.median < 0 | is.na(b.median)) & 
    (!is.na(b.egger) & b.egger < 0 | is.na(b.egger)) & 
    (!is.na(mrpresso.causal_estimate.raw) & mrpresso.causal_estimate.raw < 0 | is.na(mrpresso.causal_estimate.raw)) & 
    (!is.na(mrpresso.causal_estimate.corrected) & mrpresso.causal_estimate.corrected < 0 | is.na(mrpresso.causal_estimate.corrected)) ~ 'Pass',
  
  TRUE ~ 'Fail'
))


table(join4[c('Significance', 'Outlier_robust_estimate_test')])

# pleiotropy
join4 <- join4 %>%
  mutate(Pleiotropy_test = case_when(
    (pval.egger_intercept >= 0.05) ~ 'Pass',
    is.na(pval.egger_intercept) ~ 'NA',
    TRUE ~ 'Fail'
  ))

# count
table(join4$Pleiotropy_test)

# combine causal estimate, Heterogeneity_test, and Pleiotropy_test
join4 %>% filter(Significance == "FDR") %>% filter(Heterogeneity_test == 'Pass' | Heterogeneity_test == 'NA') %>% nrow() 
join4 %>% filter(Significance == "FDR") %>% filter(Heterogeneity_test == 'Pass' | Heterogeneity_test == 'NA')  %>% filter(Outlier_robust_estimate_test == 'Pass' | Outlier_robust_estimate_test == 'NA') 
join4 %>% filter(Significance == "FDR") %>% filter(Heterogeneity_test == 'Pass'  | Heterogeneity_test == 'NA') %>% filter(Outlier_robust_estimate_test == 'Pass' | Outlier_robust_estimate_test == 'NA') %>% filter(Pleiotropy_test == 'Pass' | Pleiotropy_test == 'NA')
join4 %>% filter(Significance == "FDR") %>% filter(Heterogeneity_test == 'Pass'  | Heterogeneity_test == 'NA') %>% filter(Outlier_robust_estimate_test == 'Pass' | Outlier_robust_estimate_test == 'NA') %>% filter(Pleiotropy_test == 'Pass' | Pleiotropy_test == 'NA') %>% filter(correct_causal_direction == TRUE) %>% nrow()
join4 %<>% mutate(Sensitivity_test = case_when(
                                         (Heterogeneity_test == 'Pass' | Heterogeneity_test == 'NA') & 
                                         (Pleiotropy_test == 'Pass'| Pleiotropy_test == 'NA') &
                                         (Outlier_robust_estimate_test == 'Pass'| Outlier_robust_estimate_test == 'NA') &
                                         correct_causal_direction == TRUE ~ 'Pass',
                                       TRUE ~ 'No'))

table(join4[c('Significance', 'Sensitivity_test')])

join4 %<>% mutate(All_pass = case_when(Significance == 'FDR' & 
                                        (Heterogeneity_test == 'FDR' | Heterogeneity_test == 'NA') & 
                                        (Pleiotropy_test == 'Pass'| Pleiotropy_test == 'NA') &
                                        (Outlier_robust_estimate_test == 'Pass'| Outlier_robust_estimate_test == 'NA') &
                                        correct_causal_direction == TRUE ~ 'Pass',
                                                 TRUE ~ 'No'))

table(join4[c('Significance', 'All_pass')])

write_tsv(join4, file = paste0('output/', cohort_name ,'.tsv'))

 
# # add name label
join4$namelabel <- NA
join4$namelabel[join4$protein %in% top10_bothsides] <- join4$shortname[join4$protein %in% top10_bothsides]

