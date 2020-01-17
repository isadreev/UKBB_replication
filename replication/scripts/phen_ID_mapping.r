#devtools::install_github("mrcieu/ieugwasr")
library(ieugwasr)
a <- gwasinfo()
b <- scan("../data/ukb-b-idlist.txt", what="character")
b <- gsub("UKB", "ukb", b)
a <- subset(a, id %in% b)
a$ukbbid <- sapply(strsplit(a$note, split=":"), function(x) x[1])

# Output the result
write.table(a, file = "../data/Phenotypes.txt",sep="\t",row.names=FALSE)
