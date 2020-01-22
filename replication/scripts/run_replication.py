# python script

import os
import argparse
import subprocess
import config

parser = argparse.ArgumentParser()
parser.add_argument('--ukbbid', help="UKB-b-id")
parser.add_argument('--dictionaryfile', help="Phenotype_Dictionary")
parser.add_argument('--phesantdir', help="")
parser.add_argument('--discoveryids', help="Discovery_IDs")
parser.add_argument('--resultdir', help="Results_Dir")
parser.add_argument('--bgen', help="Bgen")
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


# python wrapper to run this
# use subprocess module, or something similar
#Rscript make_phenotyoe_file.r $ukbbid $dictionaryfile $phesantdir $discoveryids $workdir $outdir

command = 'Rscript'
args1 = [args.ukbbid,args.dictionaryfile,args.phesantdir,args.discoveryids,'../data/']
path2script = 'make_phenotype_file.r'

#subprocess.call([command, args1, path2script], shell=True)
subprocess.Popen([command, path2script]+args1)

exit()

# this produces /path/to/results/$ukbid/phen.txt

# now run bolt on discovery data

bgenMinMaf='0.001'

numThreads=config.numThreads


#check covariate data
covarFile=row[string.uppercase.index(config.covarFile)]
if covarFile == "" and jobRun == True:
	logging.info("Missing covar file")
	sheet_connector.update(uid=job_id, column=config.statusCol, value="Missing covar file")
	mail.run(to=mailTo,subject="GWAS pipeline error - "+job_id,message="The GWAS job "+job_id+" did not complete, there was an error - covariate file name not provided")
	jobRun=False

#get covariate column info
	cCovarCols=row[string.uppercase.index(config.cCovarFileCol)]
	qCovarCols=row[string.uppercase.index(config.qCovarFileCol)]
	covarString = ''
	if jobType == 'BOLT-LMM':
		if cCovarCols != '':
			cc = cCovarCols.split(",")
			for c in cc:
				covarString += ' --covarCol='+c.strip()
		if qCovarCols != '':
			qc = qCovarCols.split(",")
			for c in qc:
				covarString += ' --qCovarCol='+c.strip()
	else:
		covList=[]
		if cCovarCols != '':
			cc = cCovarCols.split(",")
			for c in cc:
				covList.append(c.strip())
		if qCovarCols != '':
			qc = qCovarCols.split(",")
			for c in qc:
				covList.append(c.strip())
		covarString =  " --covar-name " + ','.join(str(v) for v in covList)

	if covarString == '' and jobRun == True:
		logging.info("Missing covar column data")
		sheet_connector.update(uid=job_id, column=config.statusCol, value="Missing covar column data")
		mail.run(to=mailTo,subject="GWAS pipeline error - "+job_id,message="The GWAS job "+job_id+" did not complete, there was an error - no covariate columnd data provided")
		jobRun=False

#check if covariate file exists
	com = "ssh "+sshCon+" ls "+filesDir+ "/"+covarFile
	#print com
	FNULL = open(os.devnull, 'w')
	s = subprocess.call([com],shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
	if s != 0 and jobRun == True:
		logging.info("Covariate file not found in "+filesDir)
		sheet_connector.update(uid=job_id, column=config.statusCol, value="Covariate file not found")
		mail.run(to=mailTo,subject="GWAS pipeline error - "+job_id,message="The GWAS job "+job_id+" did not complete, there was an error - missing covariate file")
		jobRun=False

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
run_bolt(bgenMinMaf,phenoFile?,args.ukbbid,1,args.resultdir,numThreads,args.phesantdir,covarFile,covarString)
# once for repl
run_bolt(bgenMinMaf,phenoFile,args.ukbbid,2,resultdir,numThreads,args.phesantdir,covarFile,covarString)

