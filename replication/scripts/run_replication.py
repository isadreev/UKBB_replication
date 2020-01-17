# python script

import argparse
import subprocess

# updazte these inputs to use argparse library
ukbbid=$1
dictionaryfile=$2
phesantdir=$3
discoveryids=$4
resultdir=$5
bgen=$6
samplefile=$7
bfile=$8


outdir="${resultdir}/${ukbbid}"
mkdir -p $outdir


# python wrapper to run this
# use subprocess module., or something similar
Rscript make_phenotyoe_file.r $ukbbid $dictionaryfile $phesantdir $discoveryids $workdir $outdir

# this produces /path/to/results/$ukbid/phen.txt

# now run bolt on discovery data


# modify this so that it knows where config is
def run_bolt(bgenMinMaf,phenoFile,phenoName,phenoCol,jobDir,numThreads,filesDir,covarFile,covarString):
	#put it together
	#" --exclude "+config.exclude + \
	#" --remove "+config.remove1 +"\\" + "\n"\
	#" --remove "+config.remove2 +"\\" + "\n"\
	com = config.bolt_exe_dir + "/bolt" +"\\" + "\n"\
	" --bfile=" + config.bfile +"\\" + "\n"\
	" --bgenFile=" + config.bgenFile + "data.chr0{1..9}.bgen" +"\\" + "\n"\
	" --bgenFile=" + config.bgenFile + "data.chr{10..22}.bgen" +"\\" + "\n"\
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
run_bolt(bgenMinMaf,phenoFile,phenoName,1,jobDir,numThreads,filesDir,covarFile,covarString)
# once for repl
run_bolt(bgenMinMaf,phenoFile,phenoName,2,jobDir,numThreads,filesDir,covarFile,covarString)

