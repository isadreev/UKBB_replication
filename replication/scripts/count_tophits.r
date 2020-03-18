library(data.table)

# Define a phenotype id
id <- "ukb-b-19953"

# This is the clump file - the independent tophits from the full 450k analysis
clumpfile <- paste0("/mnt/storage/private/mrcieu/research/scratch/IGD/data/public/", id, "/clump.txt")

# This is the discovery results - 22k SNPs tested in half of the dataset
discoveryfile <- paste0("/mnt/storage/private/mrcieu/research/mr-eve/UKBB_replication/replication/results/", id, "/discovery.statsfile.txt.gz")

# This is the replication results - 22k SNPs tested in the other half of the dataset
replicationfile <- paste0("/mnt/storage/private/mrcieu/research/mr-eve/UKBB_replication/replication/results/", id, "/replication.statsfile.txt.gz")

# This is the list of SNPs that were found to be significant in 450k samples
clump <- scan(clumpfile, what=character())

# Read in our discovery and replication data
replication <- fread(replicationfile, header=TRUE)
discovery <- fread(discoveryfile, header=TRUE)

# Find our list of SNPs in the discovery data
subset(discovery, SNP %in% clump)

# How many have p < 5e-8?
sum(discovery$P_BOLT_LMM_INF < 5e-8)

# Compare to how many in the full 450k analysis
length(clump)



# Do this for all the IDs
# Save the counts in a table
