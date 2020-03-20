#devtools::install_github("mrcieu/ieugwasr")
library(ieugwasr)
a <- gwasinfo()
b <- scan("ukb-b-idlist.txt", what="character")

# Take only BMI and CHD
#b <- c("ukb-b-19953",")
b <- b[1:5]

a <- subset(a, id %in% b)
a$ukbbid <- sapply(strsplit(a$note, split=":"), function(x) x[1])

# Output the result
write.table(a, file = "Phenotypes.txt",sep="\t",row.names=FALSE)
