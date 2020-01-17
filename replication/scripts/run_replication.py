# python script

import os
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--ukbbid', help="UKB-b-id")
parser.add_argument('--dictionaryfile', help="")
parser.add_argument('--phesantdir', help="")
parser.add_argument('--discoveryids', help="")
parser.add_argument('--resultdir', help="")
parser.add_argument('--bgen', help="")
parser.add_argument('--samplefile', help="")
parser.add_argument('--bfile', help="")
args = parser.parse_args()

print(args)


# update these inputs to use argparse library
#ukbbid=$1
#dictionaryfile=$2
#phesantdir=$3
#discoveryids=$4
#resultdir=$5
#bgen=$6
#samplefile=$7
#bfile=$8


outdir=args.resultdir + '/' + args.ukbbid
print(outdir)
os.makedirs(outdir, exist_ok=True)
exit()

# python wrapper to run this
# use subprocess module, or something similar
#Rscript make_phenotyoe_file.r $ukbbid $dictionaryfile $phesantdir $discoveryids $workdir $outdir

command = 'Rscript'
path2script = 'make_phenotype_file.r'

subprocess.call([command, args, path2script], shell=True)


# this produces /path/to/results/$ukbid/phen.txt

# now run bolt on discovery data


#paths
bioPath='/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/'

#fixed variables
LDscoresFile=bioPath+'scripts/software/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz'
geneticMapFile=bioPath+'scripts/software/BOLT-LMM_v2.3.2/tables/genetic_map_hg19_withX.txt.gz'

#static
modelSnps=bioPath+'data/model_snps_for_grm/grm6_snps.prune.in'

bolt_exe_dir=bioPath+'scripts/software/BOLT-LMM_v2.3.2'
boltSampleFile='/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample'
bfile=bioPath+'data/bolt_bfile/grm6_european_filtered_ieu'
bgenFile=bioPath+'data/dosage_bgen_filtered/'
#'/mnt/storage/private/mrcieu/research/scratch/UKBIOBANK_GWAS_Pipeline/data/scratch/snp-filtering-list-all'
bgenFile='/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/'
#boltExclude='/mnt/storage/private/mrcieu/research/scratch/UKBIOBANK_GWAS_Pipeline/data/scratch/snp-filtering-list-all/'

# modify this so that it knows where config is
def run_bolt(bgenMinMaf,phenoFile,phenoName,phenoCol,jobDir,numThreads,filesDir,covarFile,covarString):
	#put it together
	#" --exclude "+config.exclude + \
	#" --remove "+config.remove1 +"\\" + "\n"\
	#" --remove "+config.remove2 +"\\" + "\n"\
	com = config.bolt_exe_dir + "/bolt" +"\\" + "\n"\
	" --bfile=" + config.bfile +"\\" + "\n"\
	" --bgenFile=" + config.bgenFile + "data.chr0{1:9}.bgen" +"\\" + "\n"\
	" --bgenFile=" + config.bgenFile + "data.chr{10:22}.bgen" +"\\" + "\n"\
	" --bgenFile=" + config.bgenFile + "data.chrX.bgen" +"\\" + "\n"\
	" --sampleFile=" + config.boltSampleFile +"\\" + "\n"\
	" --geneticMapFile=" + config.geneticMapFile +"\\" + "\n"\
	" --bgenMinMAF=" + bgenMinMaf +"\\" + "\n"\
	" --phenoFile="+filesDir+"/" + phenoFile +"\\" + "\n"\
	" --phenoCol=" + phenoCol +"\\" + "\n"\
	" --covarFile="+filesDir+"/" + covarFile +"\\" + "\n"\
	+ covarString +"\\" + "\n"\
	" --lmm" +"\\" + "\n"\
	" --LDscoresFile=" + config.LDscoresFile +"\\" + "\n"\
	" --LDscoresMatchBp" +"\\" + "\n"\
	" --numThreads="+str(numThreads) +"\\" + "\n"\
	" --verboseStats" +"\\" + "\n"\
	" --modelSnps "+config.modelSnps +"\\" + "\n"\
	" --statsFileBgenSnps=" + jobDir +"/"+phenoName+"_imputed.txt.gz" +"\\" + "\n"\
	" --covarMaxLevels=30" +"\\" + "\n"\
	" --statsFile=" + jobDir +"/"+phenoName+"_out.txt.gz"
	return com


# Run this twice
# once for discovery
run_bolt(bgenMinMaf,phenoFile,ukbbid,1,resultdir,numThreads,phesantdir,covarFile,covarString)
# once for repl
run_bolt(bgenMinMaf,phenoFile,ukbbid,2,resultdir,numThreads,phesantdir,covarFile,covarString)

