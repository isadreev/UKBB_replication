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
../data/discoveryids.txt
```

3. Create a script which creates a phenotype for Bolt LMM

Input: 
- ukb-b id
- phenotype dictionary (step 1)
Output:
- Phenotype file ready for GWAS analysis

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
python run_replication.py \
	--ukbbid ukb-b-17314 \
	--dictionaryfile ../data/dict.rdata \
	--phesantdir /mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/derived/phesant_mod \
	--discoveryids ../data/discoveryids.txt \
	--resultdir ../results \
	--bolt_exe_dir /mnt/storage/home/kf19639/downloads/BOLT-LMM_v2.3.4 \
	--bfile /mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/bolt_bfile/grm6_european_filtered_ieu \
	--bgenDir /mnt/storage/home/kf19639/repo/UKBB_replication/replication/data \
	--boltSampleFile /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample \
	--geneticMapFile /mnt/storage/home/kf19639/downloads/BOLT-LMM_v2.3.4/tables/genetic_map_hg19_withX.txt.gz\
	--bgenMinMaf 0.001 \
	--covarFile /mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/covariates/data.covariates_ieu.bolt.txt \
	--qcovarCol '{sex,chip}' \
	--LDscoresFile /mnt/storage/home/kf19639/downloads/BOLT-LMM_v2.3.4/tables/LDSCORE.1000G_EUR.tab.gz \
	--numThreads 3 \
	--modelSnps /mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/model_snps_for_grm/grm6_snps.prune.in \
	--linker /mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/linker_app15825.csv
```


from types import SimpleNamespace 

args = SimpleNamespace(ukbbid='ukb-b-17314', dictionaryfile='../data/dict.rdata', phesantdir='/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/derived/phesant_mod', discoveryids='../data/discoveryids.txt', resultdir='../results', bolt_exe_dir= '/mnt/storage/home/kf19639/downloads/BOLT-LMM_v2.3.4', bfile='/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/bolt_bfile/grm6_european_filtered_ieu', bgenDir='/mnt/storage/home/kf19639/repo/UKBB_replication/replication/data', boltSampleFile='/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample', geneticMapFile='/mnt/storage/home/kf19639/downloads/BOLT-LMM_v2.3.4/tables/genetic_map_hg19_withX.txt.gz', bgenMinMaf='0.001', covarFile='/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/covariates/data.covariates_ieu.bolt.txt', qcovarCol='{sex,chip}', LDscoresFile='/mnt/storage/home/kf19639/downloads/BOLT-LMM_v2.3.4/tables/LDSCORE.1000G_EUR.tab.gz', 
umThreads='3', modelSnps='/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/model_snps_for_grm/grm6_snps.prune.in', linker='/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/linker_app15825.csv')

TODO:
create run_replications.py and test it works with one phenotype


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



Bolt filtering SNPs

1. Get a list of all the SNPs that were previously discovered (< 50k)

less /mnt/storage/private/mrcieu/research/scratch/IGD/data/public/ukb-b*/clump.txt | sort -u > snplist.txt

2. Extract those SNPs from the bgenfiles into a new bgenfile

https://www.well.ox.ac.uk/~gav/qctool_v2/
module add apps/qctool/2.0rc4

/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/

snplist.txt

a. from each chromosome

qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr01.bgen -incl-rsids snplist.txt -og extract.chr01.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr02.bgen -incl-rsids snplist.txt -og extract.chr02.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr03.bgen -incl-rsids snplist.txt -og extract.chr03.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr04.bgen -incl-rsids snplist.txt -og extract.chr04.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr05.bgen -incl-rsids snplist.txt -og extract.chr05.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr06.bgen -incl-rsids snplist.txt -og extract.chr06.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr07.bgen -incl-rsids snplist.txt -og extract.chr07.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr08.bgen -incl-rsids snplist.txt -og extract.chr08.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr09.bgen -incl-rsids snplist.txt -og extract.chr09.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr10.bgen -incl-rsids snplist.txt -og extract.chr10.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr11.bgen -incl-rsids snplist.txt -og extract.chr11.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr12.bgen -incl-rsids snplist.txt -og extract.chr12.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr13.bgen -incl-rsids snplist.txt -og extract.chr13.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr14.bgen -incl-rsids snplist.txt -og extract.chr14.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr15.bgen -incl-rsids snplist.txt -og extract.chr15.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr16.bgen -incl-rsids snplist.txt -og extract.chr16.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr17.bgen -incl-rsids snplist.txt -og extract.chr17.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr18.bgen -incl-rsids snplist.txt -og extract.chr18.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr19.bgen -incl-rsids snplist.txt -og extract.chr19.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr20.bgen -incl-rsids snplist.txt -og extract.chr20.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr21.bgen -incl-rsids snplist.txt -og extract.chr21.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr22.bgen -incl-rsids snplist.txt -og extract.chr22.bgen
qctool -bgen-bits 8 -g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chrX.bgen -incl-rsids snplist.txt -og extract.chrX.bgen



qctool -bgen-bits 16 -g extract.chr01.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr02.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr03.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr04.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr05.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr06.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr07.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr08.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr09.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr10.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr11.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr12.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr13.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr14.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr15.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr16.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr17.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr18.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr19.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr20.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr21.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr22.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chrX.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -og joined.bgen -os joined.sample






qctool -g extract.chr01.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr02.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr03.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr04.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr05.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr06.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr07.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr08.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr09.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr10.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr11.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr12.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr13.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr14.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr15.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr16.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr17.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr18.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr19.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr20.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr21.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -g extract.chr22.bgen -s /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample -og joined.bgen -os joined.sample





3. Only use this new bgenfile in your run_replication.py




## To run gwas on bc4

Use Snakemake. This will submit the jobs to slurm. You need to have snakemake running in a 'screen' session

To enter screen:

```
screen
```

to come out of screen, while it's still running

```
ctrl+a, d
```

to go back into screen

```
screen -r
```

to list screen sessions

```
screen -list
```


So when you're in a screen session, to run the snakemake process:

```
module add languages/anaconda3/5.2.0-tflow-1.11
snakemake -prk \
-j 400 \
--cluster-config bc4-cluster.json \
--cluster "sbatch \
  --job-name={cluster.name} \
  --partition={cluster.partition} \
  --nodes={cluster.nodes} \
  --ntasks-per-node={cluster.ntask} \
  --cpus-per-task={cluster.ncpu} \
  --time={cluster.time} \
  --mem={cluster.mem} \
  --output={cluster.output}"
```





