#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)

args <- commandArgs(T)

dir <- args[1]
exp <- args[2]
out <- args[3]

# List available GWASs
ao <- available_outcomes()

# ==================== Full MR ==============================

# Get instruments for BMI
exposure_dat <- extract_instruments(exp)

# Get effects of instruments on outcome (CHD)
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=out)

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

# Read the results of GWAS
#disc_gwas = read_exposure_data("discovery.out.txt.gz")

disc_gwas = read.table(gzfile(paste(dir,exp,"discovery.out.txt.gz",sep="")),sep="\t",header=TRUE)
repl_gwas = read.table(gzfile(paste(dir,exp,"replication.out.txt.gz",sep="")),sep="\t",header=TRUE)

# Obtain significant SNPs
index <- disc_gwas$P_BOLT_LMM < 5e-8

repl_gwas_s <- repl_gwas[index,]

out_gwas <- read.table(gzfile(paste(dir,out,"replication.out.txt.gz",sep="")),sep="\t",header=TRUE)

out_gwas_s <- out_gwas[index,]

# Harmonise the exposure and outcome data
dat <- harmonise_data(repl_gwas_s, out_gwas_s)

# Perform MR
res <- mr(dat)

write.table(res, file = "MR_repl.txt", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


