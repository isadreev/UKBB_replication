library(data.table)
require(tidyverse)
library(dplyr)

args <- commandArgs(T)
igddir <- args[1]
datadir <- args[2]
resultsdir <- args[3]

# Read all phenotype names and define each phenotype id

phen_all <- read.table(paste(datadir,"/ukb-b-idlist.txt",sep=""))

# Replace all capital letters with lowercase
phen_all <- phen_all %>% mutate(V1 = tolower(V1))

# Without using the threshold for significance

instrument_counts <- lapply(phen_all[,1], function(id) {
  tophits <- scan(file.path(igddir, id, "clump.txt"), what=character())
  
  discovery <- fread(file.path(resultsdir, id, "discovery.statsfile.txt.gz")) %>%
    subset(., SNP %in% tophits)

  replication <- fread(file.path(resultsdir, id, "replication.statsfile.txt.gz")) %>%
    subset(., SNP %in% tophits)

  dat <- inner_join(discovery, replication, by="SNP") %>% mutate(id=id)

  return(dat)

}) %>% bind_rows()

save(instrument_counts, file = paste(datadir,"/instrument_list.rdata",sep=""))


# Using the threshold for significance

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
    q <- as.numeric(ncol(discovery))
    ds <- sum(discovery[,..q] < 5e-8)
    dst <- sum((discovery[,..q] < 5e-8)&(discovery$SNP %in% clump))
  } else {
    ds <- NA
    dst <- NA
  }
  
  if (!is.na(file.size(replicationfile))) {
    replication <- fread(replicationfile, header=TRUE)
    # Find our list of SNPs in the replication data
    q <- as.numeric(ncol(replication))
    rs <- sum(replication[,..q] < 5e-8)
    rst <- sum((replication[,..q] < 5e-8)&(replication$SNP %in% clump))
  } else {
    rs <- NA
    rst <- NA
  }
  
  if (!is.na(dst)&!is.na(rst)) {
    weak <- dst-rst
    } else {
      weak <- NA
    }

  temp <- data.frame(phen=id, total=nf, sign_disc=ds, sign_repl=rs, sign_disc_top=dst, sign_repl_top=rst, weak_instr=weak)
  
  out <- rbind(out,temp)
  
} 

write.table(out, file = paste(datadir,"/count_tophits.txt",sep=""), append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

