# python script

import os
import argparse
import subprocess

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
parser.add_argument('--pcFile')
parser.add_argument('--pcCovarCol')
parser.add_argument('--LDscoresFile')
parser.add_argument('--numThreads')
parser.add_argument('--modelSnps')
parser.add_argument('--linker')
parser.add_argument('--column', help="discovery or replication")

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
def bolt_command(phenoName, bolt_exe_dir, bfile, bgenDir, boltSampleFile, geneticMapFile, bgenMinMaf, phenoFile, phenoCol, covarFile, qcovarCol, LDscoresFile, numThreads, modelSnps, resultdir, outname):
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
	" --statsFile=" + outname
	return com


def lm_command(phenoName, bolt_exe_dir, bfile, bgenDir, boltSampleFile, geneticMapFile, bgenMinMaf, phenoFile, phenoCol, pcFile, pcCovarCol, LDscoresFile, numThreads, modelSnps, resultdir, outname):
	#put it together
	com = bolt_exe_dir + "/bolt" + \
	" --bfile=" + bfile + \
	" --bgenFile=" + bgenDir + "/" + "extract.chr0{1:9}.bgen" + \
	" --bgenFile=" + bgenDir + "/" + "extract.chr{10:22}.bgen" + \
	" --sampleFile=" + boltSampleFile + \
	" --bgenMinMAF=" + bgenMinMaf + \
	" --phenoFile=" + phenoFile + \
	" --phenoCol=" + phenoCol + \
	" --covarFile=" + pcFile + \
	" --qCovarCol=" + pcCovarCol + \
	" --numThreads="+str(numThreads) + \
	" --verboseStats" + \
	" --statsFileBgenSnps=" + resultdir +"/"+phenoName+"/"+phenoCol+".statsfile.txt.gz" + \
	" --covarMaxLevels=30" + \
	" --statsFile=" + outname
	return com


outname = args.resultdir + "/" + args.ukbbid + "/" + args.column + ".out.txt.gz"

cmd = bolt_command(args.ukbbid,args.bolt_exe_dir, args.bfile, args.bgenDir, args.boltSampleFile, args.geneticMapFile, args.bgenMinMaf, outdir+"/phen.txt", args.column, args.covarFile, args.qcovarCol, args.LDscoresFile, args.numThreads, args.modelSnps, args.resultdir, outname)
subprocess.run(cmd, shell=True)


# If lmm command fails then perform simple lm instead:
if os.stat(outname).st_size == 0:
	cmd = lm_command(args.ukbbid,args.bolt_exe_dir, args.bfile, args.bgenDir, args.boltSampleFile, args.geneticMapFile, args.bgenMinMaf, outdir+"/phen.txt", args.column, args.pcFile, args.pcCovarCol, args.LDscoresFile, args.numThreads, args.modelSnps, args.resultdir, outname)
	subprocess.run(cmd, shell=True)
