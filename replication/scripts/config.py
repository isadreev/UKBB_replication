#sheet ID
#spreadsheetId = '1-3Xlwni04K2vxhRRTNfFZ2Rk8dxNQaJJxhH_5MUs1cM'
#spreadsheetId = '1-r1PxZtdI4mVZhRMF_dbTmYejZE3X93V7PIWLl5qFfc'
sheetID='1_yZ4Xkj34-y_t1DhPHWXNeLQlx_nc0dqRH_GWufGKQc'

#sheet columns
runLoc = "C"
uName = "D"
email = "E"
jobType = "F"
model = "G"
phenoName = "H"
phenoDesc = "I"
bb_cols = "J"
phenoFile="M"
phenoScript="N"
phenoFileCol="O"
covarFile="P"
cCovarFileCol="Q"
qCovarFileCol="R"
#bgenMinMafCol="K"
#runLoc = "L"
jobStatus = "S"
idCol  = "T"
subCol = "U"
comCol = "V"
statusCol = "W"
debugCol = "X"

#main path
bioPath='/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/'
rdsfPath='/projects/MRC-IEU/research/data/ukbiobank/software/gwas_pipeline/dev/release_candidate/'

#directories
phenoDir = bioPath+'data/phenotypes'

#fixed variables
LDscoresFile=bioPath+'scripts/software/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz'
geneticMapFile=bioPath+'scripts/software/BOLT-LMM_v2.3.2/tables/genetic_map_hg19_withX.txt.gz'

#bolt
bolt_exe_dir=bioPath+'scripts/software/BOLT-LMM_v2.3.2'
boltSampleFile='/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample'
bfile=bioPath+'data/bolt_bfile/grm6_european_filtered_ieu'
#bgenFile=bioPath+'data/dosage_bgen_filtered/'
#'/mnt/storage/private/mrcieu/research/scratch/UKBIOBANK_GWAS_Pipeline/data/scratch/snp-filtering-list-all'
bgenFile='/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/'
#boltExclude='/mnt/storage/private/mrcieu/research/scratch/UKBIOBANK_GWAS_Pipeline/data/scratch/snp-filtering-list-all/'

#plink
plink_exe_dir=bioPath+'scripts/plink/plink_v2.00/'
plinkSampleFile='/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22_plink.sample'
plink_ind_exclude=bioPath+'data/plink_exclusions/data.pipeline_plink_exclusions.txt'
#plink_snps_extract='/mnt/storage/private/mrcieu/research/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/UKBIOBANK_HRC_only_imputation/graded_snp_filter_lists/'
#plink_covars=bioPath+'/data/scratch/plink_500K/covariates_for_plink_pipeline.txt'
#plink_keep='/mnt/storage/home/be15516/BB_GWAS/eur.ids'

#static
#exclude=bioPath+'data/bfiles/UKBioBiLallfreqSNPexclude.dat'
modelSnps=bioPath+'data/model_snps_for_grm/grm6_snps.prune.in'
#remove1='/mnt/storage/private/mrcieu/research/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/derived/ancestry/data.less_stringent_non_europeans_exclusions.txt'
#remove2='/mnt/storage/private/mrcieu/research/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/derived/standard_exclusions/meta.recommended_exclusions.txt'

#ssh
sshCon="be15516@bc4login2.acrc.bris.ac.uk"

#Run location
configLoc="Bristol4"

#email account for sending emails
sender_email="elswob@gmail.com"

#number of nodes required to submit a job
minNode=1

#number of threads
numThreads=14

#location of log file
#logLoc="/home/be15516/mounts/rdsf_bb/data/logs/"
logLoc="/home/be15516/projects/GWAS_pipeline/BiobankGWAS_P4/"
#logLoc="/Users/be15516/projects/BiobankGWAS/logs/"