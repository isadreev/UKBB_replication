library(dplyr)
args <- commandArgs(T)

ukbbid <- args[1]
dictionaryfile <- args[2]
phesantdir <- args[3]
linkerfile <- args[4]
secondlinkerfile <- args[5]
discoveryids <- args[6]
output <- args[7]

args

# Extract from the phesant file:
# FID IID phenotype

# e.g.
# ukbbid <- "ukb-b-17314"
# dictionaryfile <- "dict.rdata"
# phesantdir <- "/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/derived/phesant_mod"

load(dictionaryfile)
library(data.table)

dir.create(output, recursive=TRUE, showWarnings=FALSE)

row <- which(dict$ukbbid == ukbbid)
a <- fread(dict$phesantfile[row], header=TRUE)
b <- subset(a, select=c("FID", dict$phesantid[row]))
names(b)[2] <- "discovery"
linker <- fread(linkerfile, header=TRUE)
load(secondlinkerfile)
linker <- merge(linker, second_linker, by.x="app", by.y="phesant")


b <- merge(b, linker, by.x="FID", by.y="phesant_mod") %>%
	{tibble(FID=.$ieu, IID=.$ieu, discovery=.$discovery)}

# read in discovery ids
dids <- scan(discoveryids, what=character())

# Set NAs for discovery and replication
b$replication <- b$discovery
b$discovery[! b$IID %in% dids] <- NA
b$replication[b$IID %in% dids] <- NA
dim(b)

outfile <- file.path(output, "phen.txt")
write.table(b, file=outfile, row=F, col=T, qu=F)
