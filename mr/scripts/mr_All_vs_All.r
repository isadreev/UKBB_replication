#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)
library(data.table)
library(gwasglue)
library(gwasvcf)
library(ieugwasr)
library(genetics.binaRies)
set_bcftools()

args <- commandArgs(T)

datadir <- args[1]
resultsdir <- args[2]
vcfdir <- args[3]
out <- args[4]

# Read all phenotype names ands define each phenotype id
phen_all <- read.table(paste(datadir,"ukb-b-idlist.txt",sep="/"))

# ==================== Full MR ==============================
  
# Extract instruments for all of the exposures
exposure_dat <- extract_instruments(phen_all[,1])

# Need to extract each of those variants from every ukb-b dataset
# First define list of unique variants
snplist <- unique(exposure_dat$SNP)
  

# Get effects of instruments on outcome using vcf files

# Get the chr:pos of every SNP
snplist_info <- ieugwasr::variants_rsid(snplist)
chrpos <- paste0(snplist_info$chr, ":", snplist_info$pos)

# Lookup from one dataset
filename <- paste(vcfdir,"/",out,"/",out,".vcf.gz",sep="")

out1 <- query_gwas(filename, chrompos=chrpos)

out2 <- gwasglue::gwasvcf_to_TwoSampleMR(out1, "outcome")

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, out2)
  
# Perform full MR
res <- mr(dat)



  
# ==================== Replication ===========================
res_r=c()
mybiglist <- list()

for (exp in phen_all[,1])
{ 
  # Read the results of GWAS
  df <- paste(resultsdir,exp,"discovery.statsfile.txt.gz",sep="/")
  dsc <- read.table(df,header=TRUE)

  rf <- paste(resultsdir,exp,"replication.statsfile.txt.gz",sep="/")
  rpc <- read.table(rf,header=TRUE)

  of <- paste(resultsdir,out,"replication.statsfile.txt.gz",sep="/")
  opt <- read.table(of,header=TRUE)


  disc_gwas <- read_exposure_data(
    filename = df,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALLELE1",
    other_allele_col = "ALLELE0",
    eaf_col = "A1FREQ",
    pval_col = tail(colnames(dsc),1)
  )
  
  repl_gwas <- read_exposure_data(
    filename = rf,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALLELE1",
    other_allele_col = "ALLELE0",
    eaf_col = "A1FREQ",
    pval_col = tail(colnames(rpc),1)
  )
  
  out_gwas <- read_outcome_data(
    filename = of,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALLELE1",
    other_allele_col = "ALLELE0",
    eaf_col = "A1FREQ",
    pval_col = tail(colnames(opt),1)
  )
  
  # Filtering SNPs for their presence in the phenotype and for p-val
  
  disc_gwas_f <- subset(disc_gwas,SNP %in% exposure_dat$SNP & pval.exposure < 5e-8)
  
  repl_gwas_f <- subset(repl_gwas,SNP %in% disc_gwas_f$SNP)
  
  out_gwas_f <- subset(out_gwas,SNP %in% disc_gwas_f$SNP)


  if (!nrow(disc_gwas_f)==0) {
    # Harmonise the exposure and outcome data
    dat_r <- harmonise_data(repl_gwas_f, out_gwas_f, action=1)
  
    # Perform MR on the replication data
    res_r <- mr(dat_r)

  } else {
    res_r <- NA
  }

  tmp <- list(Full=subset(res, id.exposure==exp), Replication=res_r)
  mybiglist[[exp]] <- tmp
}

# Save all results
save(mybiglist, file = paste(resultsdir,out,"MR_vs_All.RData",sep="/"))

