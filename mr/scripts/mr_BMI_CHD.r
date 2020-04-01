#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)
library(data.table)

args <- commandArgs(T)

dir <- args[1]
exp <- args[2]
out <- args[3]

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
    filename = paste(dir,exp,"/discovery.statsfile.txt.gz",sep=""),
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
    filename = paste(dir,exp,"/replication.statsfile.txt.gz",sep=""),
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
    filename = paste(dir,out,"/replication.statsfile.txt.gz",sep=""),
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALLELE1",
    other_allele_col = "ALLELE0",
    eaf_col = "A1FREQ",
    pval_col = "P_BOLT_LMM_INF"
)

# FIltering SNPs for their presence in the phenotype and for p-val

disc_gwas_f <- subset(disc_gwas,SNP %in% exposure_dat$SNP & pval.exposure < 5e-8)

repl_gwas_f <- subset(repl_gwas,SNP %in% disc_gwas_f$SNP)

out_gwas_f <- subset(out_gwas,SNP %in% disc_gwas_f$SNP)

# Harmonise the exposure and outcome data
dat_r <- harmonise_data(repl_gwas_f, out_gwas_f)

# Perform MR on the replication data
res_r <- mr(dat_r)



# Save all results
save(res, res_r, file = "MR_BMI_CHD.RData")


