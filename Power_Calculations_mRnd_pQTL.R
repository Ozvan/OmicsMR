####Based on script Power_Calculations_mRnd.R, but adapted for pQTL files
setwd("/PROJECTS/T2DGGI/OmicsMR/Analyses_Revision1")
library(parallel)
path.pQTL = "/PROJECTS/T2DGGI/OmicsMR/IVs/harmonized_results_aug2025/"

############################Importation of Data
#Import the list of tested genes and proteins
load("../AllTestedGenes_eQTL_pQTL_Blood.RData")

##Proportion of cases in GWAS of each ancestry
GWAS.numbers <- data.frame(Cases = c(50251,88109,242283,29375,1602,16832),
                           Controls = c(103909,339395,1569734,59368,976,33767),
                           pop = c("AFR", "EAS", "EUR", "AMR", "SAF", "SAS")) #AFR is actually AFR AMR, and AMR is HIS
rownames(GWAS.numbers) <- GWAS.numbers$pop
GWAS.numbers$Total <- GWAS.numbers$Cases+GWAS.numbers$Controls


###Summary stats for pQTL: import all the files and rbind to have a single file
deCODE_EUR <- do.call(rbind, lapply(list.files(path=paste0(path.pQTL, "pQTL_deCODE_plasma_EUR/")), function(z) read.table(paste0(path.pQTL, "pQTL_deCODE_plasma_EUR/", z), header = T, as.is = T)))
Kyoto_EAS <- do.call(rbind, lapply(list.files(path=paste0(path.pQTL, "pQTL_Kyoto_plasma_EAS/")), function(z) read.table(paste0(path.pQTL, "pQTL_Kyoto_plasma_EAS/", z), header = T, as.is = T, sep = "\t")))
ARIC_AFR <- do.call(rbind, lapply(list.files(path=paste0(path.pQTL, "pQTL_ARIC_plasma_AFR/")), function(z) read.table(paste0(path.pQTL, "pQTL_ARIC_plasma_AFR/", z), header = T, as.is = T, sep = "\t")))

ALL.pQTL <- list(EUR = deCODE_EUR, AFR = ARIC_AFR, EAS = Kyoto_EAS)
#For each file, get the protein name and compute R2
ALL.pQTL <- lapply(ALL.pQTL, function(anc){
                    anc$R2 <- 2*anc$eaf.exposure*(1-anc$eaf.exposure)*anc$beta.exposure**2
                    anc$prot.name <- sapply(anc$exposure, function(z) strsplit(z, split = "\\.")[[1]][1])
                    anc} )


############################Functions
##Function to get summary stats of the IVs: now only based on harmonized file, no need to further select the IVs
## list.harmonized.file is a list of harmonized files, one for each ancestry
function.getivs.SS <- function(gene, list.harmonized.file){ #list.clump: list with the summary stats from clump file in all tissues
  SS.pQTL <- lapply(names(list.harmonized.file), function(z) subset(list.harmonized.file[[z]], prot.name == gene))
  names(SS.pQTL) <- names(list.harmonized.file)
  SS.pQTL
}

#Function to calculate statistical power, adapted from: https://shiny.cnsgenomics.com/mRnd/
power_mRnd = function(K = 0.1, OR, N, R2, alpha = 0.05){
  bMR=K*(OR/(1+K*(OR-1))-1)
  varbMR=(K*(1-K)-bMR**2)/(N*R2)
  NCP = bMR**2/varbMR
  crit <- qchisq(1 - alpha, df = 1)
  power <- 1 - pchisq(crit, df = 1, ncp = NCP)
  return(power)
}


###Calculation of MDE
MDE_mRnd = function(N, R2, alpha = 0.05, target_power = 0.8, OR_range = c(1.01, 5), K = 0.1){
  compute_power <- function(OR) {
    bMR <- K * (OR / (1 + K * (OR - 1)) - 1)
    varbMR <- (K * (1 - K) - bMR^2) / (N * R2)
    NCP <- bMR^2 / varbMR
    crit <- qchisq(1 - alpha, df = 1)
    power <- 1 - pchisq(crit, df = 1, ncp = NCP)
    return(power - target_power)  # We want this to be zero
  }
  # Check if there's a sign change over the interval, otherwise means that the study is largely underpowered
  lower <- compute_power(OR_range[1])
  upper <- compute_power(OR_range[2])

  if (lower * upper > 0) {
    warning("No root found in OR range — power too low/high across range")
    return(NA)
  }

  # Use root-finding to solve for the OR
  result <- uniroot(compute_power, interval = OR_range)
  return(result$root)
}


#Function to wrap around the different ancestries
MDE_mRnd_wrap <- function(R2bytrait, mol_trait, ancestry){
  unique(MDE_mRnd(N = subset(GWAS.numbers, pop == ancestry)$Total,
             R2 = R2bytrait[[mol_trait]][[ancestry]] ))
}



#####################Application on pQTL 
#Get the IVs of all genes  tested
allprot.tested <- mclapply(unique(unlist(Tested.eQTL.pQTL.blood$pQTL[c("deCODE_EUR", "ARIC_AFR", "Kyotot_EAS")])), function(z) function.getivs.SS(z, list.harmonized.file = ALL.pQTL), mc.cores = 20)
##Compute R2 for each ancestry
allprot.tested.R2 <- lapply(allprot.tested, function(gene){ tt <- lapply(gene, function(z) sum(z$R2)) ; names(tt) <- names(allprot.tested[[1]]) ; tt})
names(allprot.tested.R2) <- unique(unlist(Tested.eQTL.pQTL.blood$pQTL[c("deCODE_EUR", "ARIC_AFR", "Kyotot_EAS")]))

#Compute MDE in each ancestry group
MDE.allprot.allanc <- mclapply(names(allprot.tested.R2), function(gene){
  do.call(cbind, lapply(names(allprot.tested.R2[[gene]]), function(Anc) 
    MDE_mRnd_wrap(R2bytrait = allprot.tested.R2, mol_trait = gene, ancestry = Anc)))
  }, mc.cores = 20)
MDE.allprot.allanc.df <- as.data.frame(do.call(rbind, MDE.allprot.allanc))
colnames(MDE.allprot.allanc.df) <- names(allprot.tested.R2[[1]])
rownames(MDE.allprot.allanc.df) <-  names(allprot.tested.R2) ; MDE.allprot.allanc.df$mol.trait <- rownames(MDE.allprot.allanc.df)
write.table(MDE.allprot.allanc.df, "MDE_pQTL_AllGenes_Blood_AllAnc.txt", row.names = F, col.names = T, quote = F) 
