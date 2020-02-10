# python script

import os
import argparse
import subprocess
import config

parser = argparse.ArgumentParser()
parser.add_argument('--ukbbid', help="UKB-b-id")
parser.add_argument('--dictionaryfile', help="Phenotype dictionary")
parser.add_argument('--phesantdir', help="")
parser.add_argument('--discoveryids', help="Discovery IDs")
parser.add_argument('--resultdir', help="Results dir")
parser.add_argument('--bolt_exe_dir')
parser.add_argument('--bfile')
parser.add_argument('--bgenDir')
parser.add_argument('--boltSampleFile')
parser.add_argument('--geneticMapFile')
parser.add_argument('--bgenMinMaf')
parser.add_argument('--covarFile')
parser.add_argument('--qcovarCol')
parser.add_argument('--LDscoresFile')
parser.add_argument('--numThreads')
parser.add_argument('--modelSnps')
parser.add_argument('--linker')

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


with open(outdir + "/phen.txt", "r") as f:
	print("phenotype file found")


# now run bolt on discovery data

# modify this so that it knows where config is
def bolt_command(phenoName,bolt_exe_dir, bfile, bgenDir, boltSampleFile, geneticMapFile, bgenMinMaf, phenoFile, phenoCol, covarFile, qcovarCol, LDscoresFile, numThreads, modelSnps, resultdir):
	#put it together
	com = bolt_exe_dir + "/bolt" + \
	" --bfile=" + bfile + \
	" --bgenFile=" + bgenDir + "/" + "extract.chr0{1:9}.bgen" + \
	" --bgenFile=" + bgenDir + "/" + "extract.chr{10:22}.bgen" + \
	" --sampleFile=" + boltSampleFile + \
	" --geneticMapFile=" + geneticMapFile + \
	" --bgenMinMAF=" + bgenMinMaf + \
	" --phenoFile=" + phenoFile + \
	" --phenoCol=" + phenoCol + \
	" --covarFile=" + covarFile + \
	" --qCovarCol=" + qcovarCol + \
	" --lmm" + \
	" --LDscoresFile=" + LDscoresFile + \
	" --LDscoresMatchBp" + \
	" --numThreads="+str(numThreads) + \
	" --verboseStats" + \
	" --modelSnps "+modelSnps + \
	" --statsFileBgenSnps=" + resultdir +"/"+phenoName+"/"+phenoCol+".statsfile.txt.gz" + \
	" --covarMaxLevels=30" + \
	" --statsFile=" + resultdir +"/"+phenoName+"/"+phenoCol+".out.txt.gz"
	return com

# Run this twice
# once for discovery
cmd = bolt_command(args.ukbbid,args.bolt_exe_dir, args.bfile, args.bgenDir, args.boltSampleFile, args.geneticMapFile, args.bgenMinMaf, outdir+"/phen.txt", "discovery", args.covarFile, args.qcovarCol, args.LDscoresFile, args.numThreads, args.modelSnps, args.resultdir)
subprocess.run(cmd, shell=True)

# for replication
cmd = bolt_command(args.ukbbid,args.bolt_exe_dir, args.bfile, args.bgenDir, args.boltSampleFile, args.geneticMapFile, args.bgenMinMaf, outdir+"/phen.txt", "replication", args.covarFile, args.qcovarCol, args.LDscoresFile, args.numThreads, args.modelSnps, args.resultdir)
subprocess.run(cmd, shell=True)

