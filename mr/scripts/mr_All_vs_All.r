#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)
library(data.table)

args <- commandArgs(T)

dir <- args[1]
exp <- args[2]


dir_res <- paste(dir,"results/",sep="")
datadir <- paste(dir,"data/",sep="")

# Read all phenotype names ands define each phenotype id
phen_all <- read.table(paste(datadir,"ukb-b-idlist.txt",sep=""))
phen_all <- as.data.frame(phen_all[c(5,12,18,27,31),])

mybiglist <- list()

for (out in phen_all[,1])
{
  # ==================== Full MR ==============================
  
  # Get instruments for BMI
  exposure_dat <- extract_instruments(exp)
  
  # Get effects of instruments on outcome (CHD)
  outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=out)
  
  nrow(data.frame(outcome_dat))
  
  # Harmonise the exposure and outcome data
  dat <- harmonise_data(exposure_dat, outcome_dat)
  
  # Perform full MR
  res <- mr(dat)
  
  # ==================== Replication ===========================
  
  # Read the results of GWAS
  
  disc_gwas <- read_exposure_data(
    filename = paste(dir_res,exp,"/discovery.statsfile.txt.gz",sep=""),
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALLELE1",
    other_allele_col = "ALLELE0",
    eaf_col = "A1FREQ",
    pval_col = "P_BOLT_LMM_INF"
  )
  
  
  repl_gwas <- read_exposure_data(
    filename = paste(dir_res,exp,"/replication.statsfile.txt.gz",sep=""),
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALLELE1",
    other_allele_col = "ALLELE0",
    eaf_col = "A1FREQ",
    pval_col = "P_BOLT_LMM_INF"
  )
  
  out_gwas <- read_outcome_data(
    filename = paste(dir_res,out,"/replication.statsfile.txt.gz",sep=""),
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALLELE1",
    other_allele_col = "ALLELE0",
    eaf_col = "A1FREQ",
    pval_col = "P_BOLT_LMM_INF"
  )
  
  # Filtering SNPs for their presence in the phenotype and for p-val
  
  disc_gwas_f <- subset(disc_gwas,SNP %in% exposure_dat$SNP & pval.exposure < 5e-8)
  
  repl_gwas_f <- subset(repl_gwas,SNP %in% disc_gwas_f$SNP)
  
  out_gwas_f <- subset(out_gwas,SNP %in% disc_gwas_f$SNP)
  
  # Harmonise the exposure and outcome data
  dat_r <- harmonise_data(repl_gwas_f, out_gwas_f, action=1)
  
  # Perform MR on the replication data
  res_r <- mr(dat_r)
  
  tmp <- list(Full=res, Replication=res_r)
  mybiglist[[out]] <- append(mybiglist, tmp)
  
}

# Save all results
save(mybiglist, file = paste(dir_res,exp,"MR_vs_All.RData"))

