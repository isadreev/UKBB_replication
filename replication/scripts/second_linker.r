library(data.table)
library(dplyr)

args <- commandArgs(T)
phendir <- args[1]
outfile <- args[2]

p1 <- file.path(phendir, "data", "derived", "phesant", "data-cont-10-40.txt")
p2 <- file.path(phendir, "data", "derived", "phesant_mod", "data-cont-10-40.txt")
a1 <- fread(p1)
a2 <- fread(p2)


second_linker <- tibble(phesant=a1$userID, phesant_mod=a2$FID)
save(second_linker, file=outfile)
