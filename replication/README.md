# Perform replication of ukbb gwas hits

For each `ukb-b` dataset

1. Obtain tophits
2. Split the UKBiobank phenotype into discovery and replication
3. Perform the analysis in each of those datasets

## Setup

Create a `config.json` file with relevant paths:

```json
{
    "genodir": "",
    "phendir": "",
    "igddir": ""
    "pipelinedir": "",
    "datadir": "",
    "resultsdir": "",
}
```

* `genodir` = path to ukbb bgen files
* `phendir` = path to phesant formatted phenotype files
* `igddir` = path to GWAS VCF files for ukb-b batch
* `pipelinedir` = path to original ukb-b gwas pipeline dir
* `datadir` = where to store generated data
* `resultsdir` = where to store generated results

https://github.com/MRCIEU/BiobankGWAS/blob/master/scripts/config.py
- this has all the files used for the original analysis

Genotype files:
https://github.com/MRCIEU/BiobankGWAS/blob/master/scripts/config.py#L47
https://github.com/MRCIEU/BiobankGWAS/blob/master/scripts/config.py#L44

The phenotype and genotype files have different IDs. They need to be linked using a linker file in the phesant directory
﻿
The list of individuals to include based on population are here: https://github.com/MRCIEU/BiobankGWAS/blob/master/scripts/config.py#L53

## Paper

Current version of the paper can be found here:
https://uob-my.sharepoint.com/personal/kf19639_bristol_ac_uk/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fkf19639%5Fbristol%5Fac%5Fuk%2FDocuments%2FProjects%2FUKBB%5Freplication


## To run everything on bc4

Use Snakemake to orchestrate the analysis. This will submit the jobs to slurm. 
So when you're in a screen session, to run the snakemake process:

```
module add languages/anaconda3/5.2.0-tflow-1.11
source ~/.bash_profile
snakemake -prk \
-j 100 \
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

