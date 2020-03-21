#devtools::install_github("mrcieu/ieugwasr")
library(ieugwasr)
a <- gwasinfo()
b <- scan("ukb-b-idlist.txt", what="character")

# Take only BMI and CHD
b <- c("ukb-b-19953","ukb-b-1668","ukb-b-3983","ukb-b-7436","ukb-b-16606")

a <- subset(a, id %in% b)
a$ukbbid <- sapply(strsplit(a$note, split=":"), function(x) x[1])

# Output the result
write.table(a, file = "Phenotypes.txt",sep="\t",row.names=FALSE)
