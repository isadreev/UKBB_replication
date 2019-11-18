# List of all IDs to analyse
# /mnt/storage/private/mrcieu/research/scratch/IGD/data/public/UKB-b-*/

# Map the IDs to the actual phenotype name, they should be here? 
# /mnt/storage/private/mrcieu/research/UKBIOBANK_Phenotypes_App_15825/data/derived/phesant


# PHEN = [list of phenotype names, ~2500]
PHEN = ["A", "B", "C"]


rule all:
	input: 
		expand("{phen}.disc.txt", phen=PHEN),
		expand("{phen}.repl.txt", phen=PHEN)

rule create_phenotype_files:
	input:
		"masterphenotype.txt"
	output:
		"{phen}_phenotype.txt"
	shell:
		"Rscript create_phenotype.r masterphenotype.txt {wildcards.phen}_phenotype.txt"

rule perform_disc:
	input:
		"{phen}_phenotype.txt"
	output:
		"{phen}.disc.txt"
	shell:
		"./rungwas {wildcards.phen}_phenotype.txt {wildcards.phen}.disc.txt 1"

rule perform_repl:
	input:
		"{phen}_phenotype.txt"
	output:
		"{phen}.repl.txt"
	shell:
		"./rungwas {wildcards.phen}_phenotype.txt {wildcards.phen}.repl.txt 2"


