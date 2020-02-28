library(data.table)

setwd("/mnt/storage/private/mrcieu/research/mr-eve/UKBB_replication/replication/data")

args <- commandArgs(T)
gwas_samplefile <- args[1]
randomseed <- as.numeric(args[2])
output <- args[3]

a <- fread(gwas_samplefile, skip=2)
sample_ids <- a$V1

set.seed(randomseed)
index <- sample(1:length(sample_ids), round(length(sample_ids) / 2))
discovery <- sample_ids[index]
replication <- sample_ids[-index]

write.table(discovery, file=output, row=F, col=F, qu=F)
