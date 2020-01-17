library(dplyr)

a <- read.table("Phenotypes.txt", he=T, stringsAsFactors=FALSE, comment="*")
b <- read.table("temp2", sep="\t", stringsAsFactors=FALSE, comment="*")
phesant <- lapply(1:nrow(b), function(x) {
	strsplit(b$V2[x], split=" ")[[1]]
})
names(phesant) <- b$V1

l <- list()
for(i in 1:nrow(a))
{
	message(i)
	phesantid <- a$ukbbid[i]
	col <- which(sapply(phesant, function(x) phesantid %in% x))
	phesantfile <- b$V1[col]
	l[[i]] <- tibble(phesantid=phesantid, ukbbid=a$id[i], phesantfile=phesantfile, filecolumn=col)
}

dict <- bind_rows(l)

length(unique(dict$phesantid))
length(unique(dict$ukbbid))
save(dict, file="dict.rdata")
