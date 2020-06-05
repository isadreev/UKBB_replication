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

# Read all phenotype names ands define each phenotype id
phen_all <- fread(paste(datadir,"ukb-b-idlist.txt",sep="/"), header=FALSE)
phen_all <- as.data.frame(phen_all)


# Extract instruments for all of the exposures
exposure_dat <- extract_instruments(phen_all[,1])

# Need to extract each of those variants from every ukb-b dataset
# First define list of unique variants
snplist <- unique(exposure_dat$SNP)
  

# Get effects of instruments on outcome using vcf files

# Get the chr:pos of every SNP
snplist_info <- ieugwasr::variants_rsid(snplist)
chrpos <- paste0(snplist_info$chr, ":", snplist_info$pos)

# Save all results
mybiglist <- list(chr_pos=chrpos,exp_dat=exposure_dat)
save(mybiglist, file = paste(datadir,"MR_prep.RData",sep="/"))
