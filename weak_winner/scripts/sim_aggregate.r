library(dplyr)

args <- commandArgs(T)
datadir <- args[1]

fn <- list.files(paste0(datadir,"/sim_data"), full=TRUE)
l <- list()
for(i in 1:((length(fn)-1)/4))
{
	message(i)
	load(fn[i])
	l[[i]] <- out
}

param <- bind_rows(l)
save(param, file=paste0(datadir,"/sim_overlap_nr.rdata"))

