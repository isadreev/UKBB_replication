# Read the linker file
linker <- read.table("/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/linker_app15825.csv",header=TRUE)

# Read the list of filtered genotype samples (exclusions, remove non-europeans)
excl <- read.table("/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/plink_exclusions/data.pipeline_plink_exclusions.txt")

# Find the number ids that are NOT in the exclusion list
filtered <- linker[!(linker$ieu %in% excl[,1]),"app"]


#for every phenotype
#check if the samples are in the filtered list
phen <- read.table("Phenotypes.txt",header=TRUE)

for (i in 1:nrow(phen))
{
  temp1 <- phen[i,c("id","ukbbid")]
  w <- as.character(temp1$ukbbid)
  p <- system(paste("find /mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/derived/phesant_mod -type f -print0 | xargs -0 grep -l",w ,sep=" "), intern = TRUE, input = w)
  temp2 <- cbind(temp1,path=p)
}










#phen$ukbbid <- sapply(strsplit(as.character(phen$ukbbid), split="#"), function(x) x[1])

#phen <- subset(phen, ukbbid %in% filtered)

# Create the phenotype file
out <- data.frame(IID=phen$ukbbid,FID=IID,disc=phen$value,repl=disc)

# Make half NA for disc and repl
out[round(nrow(out)/2)+1:nrow(out),"disc"] <- NA
out[1:round(nrow(out)/2),"repl"] <- NA

# Output the result
write.table(out, file = "Phenotype_name.txt")


# Create the dictionary file
dict <- data.frame(phenid=phen$id,ukbbid=phen$ukbbid,file_location=)

#find . -type f -print0 | xargs -0 grep -l "4576704"

#find /mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/derived/phesant_mod -type f -print0 | xargs -0 grep -l 4576704



ww <- system("find /mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/derived/phesant_mod -type f -print0 | xargs -0 grep -l 45767", intern = TRUE)





