# Perform replication of ukbb gwas hits


1. map ukb-b IDs to phesant ids, create a `phenotype dictionary`

```
cd scripts
Rscript phen_ID_mapping.r
./get_headers.sh
Rscript make_phenotype_dictionary.r
```

2. Choose individuals to be in discovery and replication

```
Rscript choose_discovery_ids.r \
/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample \
100 \
/path/to/results/discoveryids.txt
```

3. Create a script which creates a phenotype for Bolt LMM

Input: 
- ukb-b id
- phenotype dictionary (step 1)
Output:
- Phenotype file ready for GWAS analysis

```
Rscript 
```


/path/to/results/
				 ukb-b-2000/
				 			phen.txt
				 			disc.txt
				 			repl.txt
				 ukb-b-2001/
				 			phen.txt
				 			disc.txt
				 			repl.txt



4. Run the replication

```
python run_replications.py \
	ukb-b-2000 \
	../data/dict.rdata \
	/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/derived/phesant_mod \
	/path/to/results/discoveryids.txt \
	/path/to/results \
	/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/ \
	... all other data files needed for gwas
```



TODO:
create run_replications.sh and test it works with one phenotype


The files that go into the GWAS are:


https://github.com/MRCIEU/BiobankGWAS/blob/master/scripts/config.py
- this has all the files used for the original analysis

The phenotype files are here:
/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/derived/phesant_mod

Genotype files:
https://github.com/MRCIEU/BiobankGWAS/blob/master/scripts/config.py#L47
https://github.com/MRCIEU/BiobankGWAS/blob/master/scripts/config.py#L44

The phenotype and genotype files have different IDs. They need to be linked using this file (I think - please check):
/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/linker_app15825.csv

﻿The list of individuals to include based on population are here:
/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/plink_exclusions/data.pipeline_plink_exclusions.txt 
as defined here https://github.com/MRCIEU/BiobankGWAS/blob/master/scripts/config.py#L53

When making the phenotype files

All the phenotype IDs that we want to analyse are here:
ls -1 /mnt/storage/private/mrcieu/research/scratch/IGD/data/public | grep "UKB-b-" > ukb-b-idlist.txt
UKB-b-10001
UKB-b-10002
UKB-b-10003
UKB-b-10008
UKB-b-10011
UKB-b-1004
UKB-b-10054
UKB-b-10056
UKB-b-10061
UKB-b-10066

You can map these ID names to the UKBB ID names using:

devtools::install_github("mrcieu/ieugwasr")
library(ieugwasr)
a <- gwasinfo()
b <- scan("ukb-b-idlist.txt", what="character")
a <- subset(a, id %in% b)
a$ukbbid <- sapply(strsplit(a$note, split=":"), function(x) x[1])

Next - find which phenotype file has the ukbbid in it. Perhaps start by making a dictionary (1-many) mapping of all which files have which ukbbid 

Now you can

1. Extract the phenotype
2. Update the IDs
3. Remove non-Europeans
4. Create the phenotype file which will have columns

IID FID disc repl

Where half disc and half repl are set to NA



