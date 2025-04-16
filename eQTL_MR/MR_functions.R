library(dplyr)
#What is changed: consider the situation where Steiger filtering results in no variants left
#And where only one variant left in Steiger filtering: heterogeneity cannot be computed

################################################
#------------------- Functions -----------------
################################################
get.qtl.data <- function(exp.var){
  data <- TwoSampleMR::format_data(type = "exposure",
                                   dat = as.data.frame(exp.var),
                                   snp_col = "rs_id",
                                   beta_col = "beta",
                                   se_col = "standard_error",
                                   effect_allele_col = "effect_allele",
                                   other_allele_col = "other_allele",
                                   eaf_col = "effect_allele_frequency",
                                   pval_col = "p_value",
                                   chr_col = "chromosome",
                                   pos_col = "base_pair_location",
                                   samplesize_col = "n",
                                   gene_col = "mol.trait")
  return(data)
}

get.t2d.data <- function(filename, snps){
  outcome_dat <- TwoSampleMR::read_outcome_data(filename = filename,
                                                snps = snps,
                                                chr_col = "chromosome",
                                                pos_col = "base_pair_location",
                                                sep=",",
                                                snp_col = "rs_id",
                                                beta_col = "beta",
                                                se_col = "standard_error",
                                                effect_allele_col = "effect_allele",
                                                other_allele_col = "other_allele",
                                                eaf_col = "effect_allele_frequency",
                                                pval_col = "p_value",
                                                samplesize_col = "n")
}

add_odds_ratios <- function(dt, beta="beta", se="standard_error"){
  dt$or <- exp(dt[, get(beta)])
  dt$or_lci95 <- exp(dt[, get(beta)] - 1.96 * dt[, get(se)])
  dt$or_uci95 <- exp(dt[, get(beta)] + 1.96 * dt[, get(se)])
  return(dt)
}

run.MRPresso <- function(my.data, nbdist=1000, beta.out="beta.outcome", beta.exp="beta.exposure", se.out="se.outcome", se.exp="se.exposure"){
  dt <- MRPRESSO::mr_presso(BetaOutcome=beta.out, 
                            BetaExposure=beta.exp, 
                            SdOutcome=se.out, 
                            SdExposure=se.exp, 
                            OUTLIERtest=TRUE, 
                            DISTORTIONtest=TRUE, 
                            data=my.data, 
                            NbDistribution=nbdist,  
                            SignifThreshold=0.05)
  #Results of the outlier tests
  #Also keep the disortion p-value and the number of snps
  res.presso <- data.table::as.data.table(dt$`Main MR results`)[`MR Analysis`=="Outlier-corrected"]
  res.presso <- add_odds_ratios(res.presso, "Causal Estimate", "Sd")
  
  res.presso <- res.presso[, .(b=`Causal Estimate`, se=`Sd`, or, or_lci95, or_uci95, pval=`P-value`,
                               method="MR_PRESSO_IVW_Outliers", exposure=my.data$exposure[1], outcome=my.data$outcome[1],
                               Distortion_pval = dt$`MR-PRESSO results`$`Distortion Test`$Pvalue,
                               nsnp = nrow(my.data)-length(dt$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`))]

  #Add the raw test
  res.presso.raw <- data.table::as.data.table(dt$`Main MR results`)[`MR Analysis`=="Raw"]
  res.presso.raw <- add_odds_ratios(res.presso.raw, "Causal Estimate", "Sd")

  res.presso.raw <- res.presso.raw[, .(b=`Causal Estimate`, se=`Sd`, or, or_lci95, or_uci95, pval=`P-value`,
                               method="MR_PRESSO_IVW_Raw", exposure=my.data$exposure[1], outcome=my.data$outcome[1],
                               nsnp = nrow(my.data))]

  # res.presso<- run.SensitivityAnalysis(my.data, res.presso)
  return(rbind(res.presso, res.presso.raw, fill = T))
}

run.SensitivityAnalysis <- function(data, res, plt=TRUE, het=TRUE, presso=TRUE, steiger = TRUE, Radial){
  fstat_overall <- mean((data$beta.exposure)^2/(data$se.exposure)^2)


  res[, Fstat:=fstat_overall]
  if(het){
    # check heterogeneity with Cochran's Q-statistic
    het <- data.table::as.data.table(TwoSampleMR::mr_heterogeneity(data)) 
    res <- merge(res, het[, .(exposure, outcome, method, Q, Q_df, Q_pval)],
                 by=c("exposure", "outcome", "method"), all.x=TRUE)
  }
  if(plt){
    # Assumption 2: check pleiotropy with MR-Egger intercept
    plt <- data.table::as.data.table(TwoSampleMR::mr_pleiotropy_test(data))  #MR Egger intercept for directional pleiotropy
    res[method=="MR Egger", `:=` (egger_intercept=plt$egger_intercept, egger_intercept_se=plt$se, egger_intercept_pval=plt$pval)]
  }
  if(presso){
    res.presso <- tryCatch(
      run.MRPresso(data), 
      error=function(e) e
    )
    
    if(inherits(res.presso, "error")) {
      print("MR-PRESSO failed. Probably not enough SNPs.")
    }else{
      res <- rbind(res, res.presso, fill=TRUE)
    }
  }

  if(steiger){
    #Run the Steiger filtering
    data.steiger <- steiger_filtering(data)
    data.steiger <- subset(data.steiger, steiger_dir & steiger_pval<0.05)
    if(nrow(data.steiger)>0){
      res.steiger <- TwoSampleMR::mr(data.steiger, method = c("mr_ivw", "mr_wald_ratio"))
      res.steiger <- data.table::as.data.table(TwoSampleMR::generate_odds_ratios(res.steiger))
      res.steiger[, `:=`(lo_ci=NULL, up_ci=NULL, id.exposure=NULL, id.outcome=NULL)]
      if(nrow(data.steiger)>1){
        #Compute heterogeneity
        het <- data.table::as.data.table(TwoSampleMR::mr_heterogeneity(data.steiger, method = "mr_ivw"))
        res.steiger <- merge(res.steiger, het[, .(exposure, outcome, method, Q, Q_df, Q_pval)],
                     by=c("exposure", "outcome", "method"), all.x=TRUE)
        res.steiger$method <- "Steiger Inverse variance weighted"
      }else{
        res.steiger$method <- "Steiger Wald ratio"
      }
    }else{
      res.steiger <- data.frame(method = "Steiger Inverse variance weighted")
    }
    res <- rbind(res, res.steiger, fill = TRUE)
  }

  if(Radial){
    #Run radial filtering, esp. for metabolites
    #Can have warning or error if model does not converge
    radial.outliers <- tryCatch(ivw_radial(data, alpha = 0.05, weights = 1), error = identity, warning = identity)
    #If there are SNPs left after filtering, proceed with the test
    if(is(radial.outliers, "error") | is(radial.outliers, "warning")){
      res.radial <- data.table(method = "Radial Inverse variance weighted", nsnp = NA)
    }else{
      if(nrow(data) == length(radial.outliers$outliers$SNP)){
        res.radial <- data.table(method = "Radial Inverse variance weighted", nsnp = 0)
      }else{
        if(is.null(ncol(radial.outliers$outliers))){
          #No outlier detected
          res.radial <- data.table(method = "Radial Inverse variance weighted", nsnp = nrow(data))
        }else{
          data.radial <- subset(data, !(SNP %in% radial.outliers$outliers$SNP))
          res.radial <- TwoSampleMR::mr(data.radial, method = "mr_ivw")
          res.radial <- data.table::as.data.table(TwoSampleMR::generate_odds_ratios(res.radial))
          res.radial[, `:=`(lo_ci=NULL, up_ci=NULL, id.exposure=NULL, id.outcome=NULL)]
          #Compute heterogeneity
          het <- data.table::as.data.table(TwoSampleMR::mr_heterogeneity(data.radial, method = "mr_ivw"))
          res.radial <- merge(res.radial, het[, .(exposure, outcome, method, Q, Q_df, Q_pval)],
                 by=c("exposure", "outcome", "method"), all.x=TRUE)
          res.radial$method <- "Radial Inverse variance weighted"
        }
      }
    }
    res <- rbind(res, res.radial, fill = TRUE)
  }


  return(res)
}

run.TwoSampleMR <- function(data, Radial){
  if(nrow(data)==0)  {
    print("No SNP in common between exposure and outcome :(")
    tmp.dt <- data.table::data.table()
    return(tmp.dt)
  }
  else if (nrow(data)==1) {
    fstat <- (data$beta.exposure)^2/(data$se.exposure)^2
    print("Only one SNP in common between exposure and outcome, running Wald ratio method")
    res <- TwoSampleMR::mr(data, method_list="mr_wald_ratio")
    res <- data.table::as.data.table(TwoSampleMR::generate_odds_ratios(res))
    res[, `:=`(lo_ci=NULL, up_ci=NULL, id.exposure=NULL, id.outcome=NULL)]
    
    tmp.dt <- data.table::as.data.table(res) %>% .[, `:=`(Fstat=fstat,Q=as.numeric(NA),Q_df=as.numeric(NA),Q_pval=as.numeric(NA),plt.egger_intercept=as.numeric(NA),plt.pval=as.numeric(NA),plt.se=as.numeric(NA))]
    return(tmp.dt)
  }
  else {
    print("Multiple SNPs in common between exposure and outcome, running MR")
    res <- TwoSampleMR::mr(data, method_list=c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
    res <- data.table::as.data.table(TwoSampleMR::generate_odds_ratios(res))
    res[, `:=`(lo_ci=NULL, up_ci=NULL, id.exposure=NULL, id.outcome=NULL)]
    
    # Run sensitivity analysis
    res <- run.SensitivityAnalysis(data, res, Radial = Radial)
    
    return(res)
  }
}

wrap.MR <- function(exposure.lst, t2d.filename, ivs.path, data.path, output.file, output.harmonized, Radial = F){
  dt <- data.table::data.table()
  
  #----------- Get outcome -----------
  ivs <- data.table::fread(ivs.path)
  outcome <- get.t2d.data(paste0(data.path, t2d.filename), snps=ivs$rs_id)
  outcome$outcome <- "T2D"
  
  for (exp in exposure.lst){
    exp.name <- stringr::str_extract(exp, "[^_]*__[^_]*__[^_]*__[^_]*__[^_]*")
    exp.build <- sub("\\..*", "", sub(".*_", "", exp))
    
    #----------- Read and clump exposure -----------
    exposure.all <- data.table::fread(paste0(data.path, exp))
    
    for (cur.trait in unique(exposure.all$mol.trait)){
      exposure <- get.qtl.data(exposure.all[mol.trait==cur.trait])
      exposure$exposure <- paste(sub("clumped__", "", exp.name), cur.trait, sep="__")
             
      #----------- Harmonize data -----------
      data <- tryCatch(
        TwoSampleMR::harmonise_data(exposure_dat=exposure, outcome_dat=outcome),
        error=function(e) e
      )
      
      if(inherits(data, "error")) {
        print("Harmonization failed.")
        next 
      }
      
      if(nrow(subset(data, mr_keep==TRUE))==0) {
        print("No SNP left after data harmonization, not able to run MR :(")
        next
      }
      data$id.exposure <- data$exposure
      data$id.outcome <- "T2D"      
      data <- subset(data, mr_keep==TRUE)
      
      #----------- Run TwoSampleMR with sensitivity analyses -----------
      res.mr <- run.TwoSampleMR(data, Radial = Radial)
      dt <- rbind(dt, res.mr, fill=TRUE)

      #Compute the I2 statistics
      dt$I2 <- (dt$Q-dt$Q_df)/dt$Q

      ###Save the harmonized data
      fwrite(data, paste0(output.harmonized, data$gene.exposure[1], ".txt"), na = "NA", sep = "\t", quote = F)
    }
  }
  
  ##Check if MR-Egger has been run for at least one molecular trait
  #Otherwise, columns are not present
  if("egger_intercept_se" %in% colnames(dt)){
	dt[method %in% c("Inverse variance weighted", "Wald ratio"), p.ivw.adj.fdr:=p.adjust(pval, method = "fdr")] %>% setnames(., old=c("b", "se", "pval", "or", "or_lci95", "or_uci95", "egger_intercept_se", "egger_intercept_pval", "p.ivw.adj.fdr"),
													 new=c("beta", "standard_error", "p_value", "odds_ratio", "ci_lower", "ci_upper", "egger_intercept_standard_error", "egger_intercept_pvalue", "p_value_ivw_fdr"))
  }else{
	dt[method %in% c("Inverse variance weighted", "Wald ratio"), p.ivw.adj.fdr:=p.adjust(pval, method = "fdr")] %>% setnames(., old=c("b", "se", "pval", "or", "or_lci95", "or_uci95", "p.ivw.adj.fdr"),
                                                                                                         new=c("beta", "standard_error", "p_value", "odds_ratio", "ci_lower", "ci_upper", "p_value_ivw_fdr"))
  }
  data.table::fwrite(dt, output.file, na = "NA", sep = "\t", quote = F)

}

