#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)

# List available GWASs
ao <- available_outcomes()

# ==================== Full MR ==============================

# Get instruments for BMI
exposure_dat <- extract_instruments("ukb-b-19953")

# Get effects of instruments on outcome (CHD)
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ukb-b-1668")

nrow(data.frame(outcome_dat))

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
res <- mr(dat)

write.table(res, file = "MR_full.txt", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


# ==================== Replication ===========================

# Read the results og GWAS
disc_gwas = read.table(gzfile("discovery.out.txt.gz"),sep="\t",header=TRUE)
repl_gwas = read.table(gzfile("replication.out.txt.gz"),sep="\t",header=TRUE)

index <- disc_gwas$P_BOLT_LMM < 5e-8


repl_gwas[index,]





res <- mr(dat)

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
res <- mr(dat)

write.table(res, file = "MR_repl.txt", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


