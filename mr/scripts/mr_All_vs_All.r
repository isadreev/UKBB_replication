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
  print(paste("exposure=",exp))
  print(paste("outcome=",out))

  # Read the results of GWAS
  df <- paste(resultsdir,exp,"discovery.statsfile.txt.gz",sep="/")
  dsc <- read.table(df,header=TRUE)

  dsc1 <- dsc[,c("SNP","ALLELE1","ALLELE0","A1FREQ","BETA", "SE",tail(colnames(dsc),1))]
  colnames(dsc1) <- c("SNP",
    "effect_allele.exposure",
    "other_allele.exposure",
    "eaf.exposure",
    "beta.exposure",
    "se.exposure",
    "pval.exposure")
  dsc2 <- dsc1[order(dsc1$SNP, dsc1$pval.exposure), ]
  dsc2 <- dsc2[ !duplicated(dsc2$SNP), ]  
  dsc3 <- cbind(dsc2,
    "exposure"=rep("exposure",nrow(dsc2)),
    "mr_keep.exposure"=rep("TRUE", nrow(dsc2)),
    "pval_origin.exposure"=rep("reported", nrow(dsc2)),
    "id.exposure"=rep(exp, nrow(dsc2)),
    "data_source.exposure"=rep("textfile", nrow(dsc2)))




  rf <- paste(resultsdir,exp,"replication.statsfile.txt.gz",sep="/")
  rpc <- read.table(rf,header=TRUE)

  rpc1 <- rpc[,c("SNP","ALLELE1","ALLELE0","A1FREQ","BETA", "SE",tail(colnames(rpc),1))]
  colnames(rpc1) <- c("SNP",
    "effect_allele.exposure",
    "other_allele.exposure",
    "eaf.exposure",
    "beta.exposure",
    "se.exposure",
    "pval.exposure")
  rpc2 <- rpc1[order(rpc1$SNP, rpc1$pval.exposure), ]
  rpc2 <- rpc2[ !duplicated(rpc2$SNP), ]  
  rpc3 <- cbind(rpc2,
    "exposure"=rep("exposure",nrow(rpc2)),
    "mr_keep.exposure"=rep("TRUE", nrow(rpc2)),
    "pval_origin.exposure"=rep("reported", nrow(rpc2)),
    "id.exposure"=rep(exp, nrow(rpc2)),
    "data_source.exposure"=rep("textfile", nrow(rpc2)))



  of <- paste(resultsdir,out,"replication.statsfile.txt.gz",sep="/")
  opt <- read.table(of,header=TRUE)

  opt1 <- opt[,c("SNP","ALLELE1","ALLELE0","A1FREQ","BETA", "SE",tail(colnames(opt),1))]
  colnames(opt1) <- c("SNP",
    "effect_allele.outcome",
    "other_allele.outcome",
    "eaf.outcome",
    "beta.outcome",
    "se.outcome",
    "pval.outcome")
  opt2 <- opt1[order(opt1$SNP, opt1$pval.outcome), ]
  opt2 <- opt2[ !duplicated(opt2$SNP), ]  
  opt3 <- cbind(opt2,
    "outcome"=rep("outcome",nrow(opt2)),
    "mr_keep.outcome"=rep("TRUE", nrow(opt2)),
    "pval_origin.outcome"=rep("reported", nrow(opt2)),
    "id.outcome"=rep(out, nrow(opt2)),
    "data_source.outcome"=rep("textfile", nrow(opt2)))


  disc_gwas <- dsc3
  repl_gwas <- rpc3
  out_gwas <- opt3



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
  
  #repl_gwas <- read_exposure_data(
  #  filename = rf,
  #  sep = "\t",
  #  snp_col = "SNP",
  #  beta_col = "BETA",
  #  se_col = "SE",
  #  effect_allele_col = "ALLELE1",
  #  other_allele_col = "ALLELE0",
  #  eaf_col = "A1FREQ",
  #  pval_col = tail(colnames(rpc),1)
  #)
  
  #out_gwas <- read_outcome_data(
  #  filename = of,
  #  sep = "\t",
  # snp_col = "SNP",
  #  beta_col = "BETA",
  #  se_col = "SE",
  #  effect_allele_col = "ALLELE1",
  #  other_allele_col = "ALLELE0",
  #  eaf_col = "A1FREQ",
  #  pval_col = tail(colnames(opt),1)
  #)
  
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

  #tmp <- list(Full=subset(res, id.exposure==exp), Replication=res_r)
  #mybiglist[[exp]] <- tmp
}

# Save all results
save(mybiglist, file = paste(resultsdir,out,"MR_vs_All.RData",sep="/"))

