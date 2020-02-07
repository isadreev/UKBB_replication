
# define the list of datasets that we want to make

with open("../data/ukb-b-idlist.txt", 'r') as f:
	DATA = f.readlines()

DATA = [x.strip().lower() for x in DATA]

DATA=DATA[0:3]


rule all:
	input: 
		expand("../results/{data}/discovery.statsfile.txt.gz", data=DATA),
		expand("../results/{data}/replication.statsfile.txt.gz", data=DATA)


rule run_gwas:
	output:
		"../results/{data}/discovery.statsfile.txt.gz",
		"../results/{data}/replication.statsfile.txt.gz"
	shell:
"""
python run_replication.py \
	--ukbbid {data} \
	--dictionaryfile ../data/dict.rdata \
	--phesantdir /mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/derived/phesant_mod \
	--discoveryids ../data/discoveryids.txt \
	--resultdir ../results \
	--bolt_exe_dir /mnt/storage/home/kf19639/downloads/BOLT-LMM_v2.3.4 \
	--bfile /mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/bolt_bfile/grm6_european_filtered_ieu \
	--bgenDir /mnt/storage/home/kf19639/repo/UKBB_replication/replication/data \
	--boltSampleFile /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample \
	--geneticMapFile /mnt/storage/home/kf19639/downloads/BOLT-LMM_v2.3.4/tables/genetic_map_hg19_withX.txt.gz \
	--bgenMinMaf '0.001' \
	--covarFile /mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/covariates/data.covariates_ieu.bolt.txt \
	--qcovarCol '{sex,chip}' \
	--LDscoresFile /mnt/storage/home/kf19639/downloads/BOLT-LMM_v2.3.4/tables/LDSCORE.1000G_EUR.tab.gz \
	--numThreads 3 \
	--modelSnps /mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/model_snps_for_grm/grm6_snps.prune.in \
	--linker /mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/linker_app15825.csv
"""

