args <- commandArgs(T)

ukbbid <- args[1]
dictionaryfile <- args[2]
phesantdir <- args[3]
discoveryids <- args[4]
output <- args[5]

# Extract from the phesant file:
# FID IID phenotype

# e.g.
# ukbbid <- "ukb-b-17314"
# dictionaryfile <- "dict.rdata"
# phesantdir <- "/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/derived/phesant_mod"

load(dictionaryfile)
library(data.table)

row <- which(dict$ukbbid == ukbbid)
a <- fread(dict$phesantfile[row], header=TRUE)
b <- subset(a, select=c("FID", "IID", dict$phesantid[row]))


# read in discovery ids
dids <- scan(discoveryids, what=character())

# Set NAs for discovery and replication
names(b)[3] <- "discovery"
b$replication <- b$discovery
b$discovery[! b$IID %in% dids] <- NA
b$replication[b$IID %in% dids] <- NA

# Here we assume all discovery are 1, however there are 2854 elements (0.6%) where it is NA or 2:
#c=as.data.frame(b)
#c[!c[,3]==1,]

outfile <- file.path(output, "phen.txt")
write.table(b, file=outfile, row=F, col=F, qu=F)
