library(data.table)
require(tidyverse)

args <- commandArgs(T)
igddir <- args[1]
datadir <- args[2]
resultsdir <- args[3]

# Read all phenotype names ands define each phenotype id
#id <- "ukb-b-19953"

phen_all <- read.table(paste(datadir,"/ukb-b-idlist.txt",sep=""))

# Replace all capital letters with lowercase
phen_all <- phen_all %>% mutate(V1 = tolower(V1))

out <- c()

for (id in phen_all[,1])
{
  # This is the clump file - the independent tophits from the full 450k analysis
  clumpfile <- paste0(igddir, "/", id, "/clump.txt")
  
  # This is the discovery results - 22k SNPs tested in half of the dataset
  discoveryfile <- paste0(resultsdir, "/", id, "/discovery.statsfile.txt.gz")
  
  # This is the replication results - 22k SNPs tested in the other half of the dataset
  replicationfile <- paste0(resultsdir, "/", id, "/replication.statsfile.txt.gz")
  
  
  # This is the list of SNPs that were found to be significant in 450k samples
  if (!is.na(file.size(clumpfile))) {
    clump <- scan(clumpfile, what=character())
    nf <- length(clump)
  } else {
    nf <- NA
  }
  
  # Read in our discovery and replication data
  if (!is.na(file.size(discoveryfile))) {
    discovery <- fread(discoveryfile, header=TRUE)
    # Find our list of SNPs in the discovery data
    subset(discovery, SNP %in% clump)
    nd <- sum(discovery$P_BOLT_LMM_INF < 5e-8)
  } else {
    nd <- NA
  }
  
  if (!is.na(file.size(replicationfile))) {
    replication <- fread(replicationfile, header=TRUE)
  }
  
  temp <- data.frame(phen=id, full=nf, replication=nd)
  
  out <- rbind(out,temp)
  
} 

write.table(out, file = paste(datadir,"/count_tophits.txt",sep=""), append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

