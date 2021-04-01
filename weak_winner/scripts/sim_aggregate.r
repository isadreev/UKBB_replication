library(dplyr)

args <- commandArgs(T)

fn <- list.files("../data/sim_data", full=TRUE)
l <- list()
for(i in 1:((length(fn)-1)/4))
{
	message(i)
	load(fn[i])
	l[[i]] <- out
}

param <- bind_rows(l)
save(param, file="../results/sim_overlap_nr.rdata")

