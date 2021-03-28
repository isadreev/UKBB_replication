library(simulateGP)
library(tidyverse)
library(parallel)

runsim <- function(p)
{
	map <- tibble(snp=1:p$nsnp, af=runif(p$nsnp, 0.01, 0.99))
	gp <- generate_gwas_params(map, h2=p$h2, S = p$S) %>%
		mutate(
			rsq=2 * beta^2 * af * (1-af)
		)
	ss <- generate_gwas_ss(gp, nid=p$nid) %>%
		mutate(tfval = (gp$beta/se)^2)
	p$nsig <- sum(ss$pval < p$pval)
	p$nsig_weak <- sum(ss$pval < p$pval & ss$tfval < 20)
	p$nsig_overestimate <- sum(ss$pval < p$pval & (abs(ss$bhat) - 1.96 * ss$se) > abs(gp$beta))
	p$min_rsq <- subset(ss, pval < p$pval) %>% {min(.$rsq)}
	return(p)
}

param <- expand.grid(
	nid=seq(10000, 250000, by=10000),
	nsnp=c(500, 1000, 5000, 10000, 50000, 100000),
	h2=c(0.1, 0.5, 0.9),
	pval=5e-8,
	S=c(-1, 0, 1),
	rep=1:10
) %>% as_tibble()
str(param)

l <- mclapply(1:nrow(param), function(i) {
	runsim(param[i,])
	}, mc.cores=16) %>%
	bind_rows()

save(l, file="../results/wc_sim.rdata")
