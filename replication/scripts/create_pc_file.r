library(dplyr)
library(data.table)

args <- commandArgs(T)

pcfile <- args[1]
exclusionfile <- args[2]
covfile <- args[3]
outfile <- args[4]

pcs <- fread(pcfile)
exclusions <- fread(exclusionfile, header=FALSE)
covs <- fread(covfile, header=TRUE)

pcs <- subset(pcs, !V1 %in% exclusions$V1)
pcs <- subset(pcs, select=-c(V2))
covs <- merge(covs, pcs, by.x="FID", by.y="V1")
names(covs)[3:44] <- paste0("V", 1:42)
write.table(covs, file=outfile, row=F, col=F, qu=F)
