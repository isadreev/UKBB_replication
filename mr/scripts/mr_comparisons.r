library(tidyverse)
library(jsonlite)
library(data.table)

config <- read_json("config.json")
datadir <- config['datadir']
resultsdir <- config['resultsdir']


phen_all <- scan(file.path(datadir,"ukb-b-idlist.txt"), what=character())

mr_res <- list()
het_res <- list()

for (id in phen_all)
{
	message(id)
	# Reading the results for each trait
 	mr_file_all <- file.path(resultsdir, id, "MR_All_vs_All_all.txt")
 	mr_file_weak <- file.path(resultsdir, id, "MR_All_vs_All_weak.txt")
 	het_file_all <- file.path(resultsdir, id, "MR_Het_All_vs_All_all.txt")
 	het_file_weak <- file.path(resultsdir, id, "MR_Het_All_vs_All_weak.txt")

	mr_table_all <- fread(mr_file_all, header=TRUE)
	mr_table_weak <- fread(mr_file_weak, header=TRUE)
	het_table_all <- fread(het_file_all, header=TRUE)
	het_table_weak <- fread(het_file_weak, header=TRUE)

	mr_table_all <- mr_table_all[complete.cases(mr_table_all[ , "b"]), ] %>% mutate(inst = "all")
	mr_table_weak <- mr_table_weak[complete.cases(mr_table_weak[ , "b"]), ] %>% mutate(inst = "weak")
	het_table_all <- het_table_all[complete.cases(het_table_all[ , "Q"]), ] %>% mutate(inst = "all")
	het_table_weak <- het_table_weak[complete.cases(het_table_weak[ , "Q"]), ] %>% mutate(inst = "weak")

	mr_res[[id]] <- bind_rows(mr_table_all, mr_table_weak)
	het_res[[id]] <- bind_rows(het_table_all, het_table_weak)
}

mr_res <- bind_rows(mr_res)
het_res <- bind_rows(het_res)

save(mr_res, het_res, file="mr_results.rdata")

mr_res %>%
	filter(method %in% c("Inverse varianse weighted", "Wald ratio")) %>%
	group_by(data, inst) %>%
	dplyr::select(id.exposure, id.outcome, b, dir) %>%
	pivot_wider(names_from=dir, values_from=b) %>%
	summarise(r = cor(AB, BA, use="pair"))

mr_res %>%
	filter(method %in% c("Inverse varianse weighted", "Wald ratio")) %>%
	mutate(s=sign(b)) %>%
	group_by(data, inst) %>%
	dplyr::select(id.exposure, id.outcome, s, dir) %>%
	pivot_wider(names_from=dir, values_from=s) %>%
	summarise(r = cor(AB, BA, use="pair"))


mr_res %>%
	filter(method %in% c("Inverse varianse weighted", "Wald ratio") & dir=="AB") %>%
	group_by(inst) %>%
	dplyr::select(id.exposure, id.outcome, pval, data) %>%
	pivot_wider(names_from=data, values_from=pval) %>%
	summarise(
		rdsig = sum(RD < 0.05, na.rm=T) / sum(!is.na(RD)),
		rsig = sum(R < 0.05, na.rm=T) / sum(!is.na(R)),
		drsig = sum(DR < 0.05, na.rm=T) / sum(!is.na(DR)),
		dsig = sum(D < 0.05, na.rm=T) / sum(!is.na(D)),
		rcompfdr = sum(DR > 0.05 & R < 0.05, na.rm=T) / sum(DR > 0.05, na.rm=T),
		dcompfdr = sum(DR > 0.05 & D < 0.05, na.rm=T) / sum(DR > 0.05, na.rm=T),
		rdcompfdr = sum(DR > 0.05 & RD < 0.05, na.rm=T) / sum(DR > 0.05, na.rm=T),
		rcompp = sum(DR < 0.05 & R > 0.05, na.rm=T) / sum(DR < 0.05, na.rm=T),
		dcompp = sum(DR < 0.05 & D > 0.05, na.rm=T) / sum(DR < 0.05, na.rm=T),
		rdcompp = sum(DR < 0.05 & RD > 0.05, na.rm=T) / sum(DR < 0.05, na.rm=T)
	) %>%
	pivot_longer(!inst, names_to="data", values_to="estimate")


mr_res %>%
	filter(method %in% c("Inverse varianse weighted", "Wald ratio") & dir=="AB") %>%
	mutate(z=b/se) %>%
	group_by(inst) %>%
	dplyr::select(id.exposure, id.outcome, z, data) %>%
	pivot_wider(names_from=data, values_from=z) %>%
	summarise(
		rcomp = cor(R, RD, use="pair"),
		dcomp = cor(R, D, use="pair"),
		rdcomp = cor(R, DR, use="pair")
	) %>%
	pivot_longer(!inst, names_to="data", values_to="estimate")

