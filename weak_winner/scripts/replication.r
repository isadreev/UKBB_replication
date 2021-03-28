library(simulateGP)
library(TwoSampleMR)
library(dplyr)
library(ggplot2)
library(parallel)
library(knitr)
args <- commandArgs(T)
datadir <- args[1]


## Outline of simulation

## 1. Create x and y in a population
## 2. Find instruments for x in discovery and test for significance
## 3. Retest in replication
## 4. Extract those effects from outcome GWAS, which could be independent, discovery or replication sample


### Run MR simulations for different sampling strategies
umvcue <- function(d_beta, r_beta, d_se, r_se, pthresh)
{
  t <- qnorm(1 - pthresh/2)
  MLE <- (r_se^2 * d_beta + d_se^2 * r_beta) / (r_se^2 + d_se^2)
  BB <- (sqrt(r_se^2 + d_se^2) / d_se^2) * (MLE - t * d_se)
  b <- MLE - (r_se^2 / sqrt(r_se^2 + d_se^2) * dnorm(BB) / pnorm(BB))
  return(b)
}

gwas_bias <- function(gwas)
{
  # sum((gwas$bhat - gwas$b)^2)
  mean(abs(gwas$bhat - gwas$b))
}
gwas_sim <- function(nid = 9000, nsnp = 20, bgx = 0.08, bxy = 0.2, buy = 0.4, bux = 0.4, out_index = 1:3000)
{
  # Simulate population
  g <- make_geno(nid, nsnp, 0.5)
  u <- rnorm(nid)
  geneff <- rep(bgx, nsnp)
  x <- make_phen(c(geneff, bux), cbind(g, u))
  y <- make_phen(c(buy, bxy), cbind(u, x))
  # Get observational estimate
  obs_beta <- summary(lm(y ~ x))$coefficients[2,1]
  # GWAS in samples
  disc_index <- 1:3000
  rep_index <- 3001:6000
  disc_gwas <- gwas(x[disc_index], g[disc_index,])
  rep_gwas <- gwas(x[rep_index], g[rep_index,])
  out_gwas <- gwas(y[out_index], g[out_index,])
  disc_gwas$b <- geneff
  rep_gwas$b <- geneff
  out_gwas$b <- geneff * bxy
  bias <- expand.grid(sample=c("disc", "rep", "out"), what=c("all", "sig"), bias=NA)
  bias$sse[1] <- gwas_bias(disc_gwas)
  bias$sse[2] <- gwas_bias(rep_gwas)
  bias$sse[3] <- gwas_bias(out_gwas)
  # what is significant
  index <- disc_gwas$pval < 5e-8
  bias$sse[4] <- gwas_bias(disc_gwas[index,])
  bias$sse[5] <- gwas_bias(rep_gwas[index,])
  bias$sse[6] <- gwas_bias(out_gwas[index,])
  bias$nsig <- sum(index)
  # Perform MR
  d <- list()
  d$disc_all <- simulateGP::merge_exp_out(disc_gwas, out_gwas)
  d$rep_all <- simulateGP::merge_exp_out(rep_gwas, out_gwas)
  d$disc_sig <- simulateGP::merge_exp_out(disc_gwas[index,], out_gwas[index,])
  d$rep_sig <- simulateGP::merge_exp_out(rep_gwas[index,], out_gwas[index,])
  d$umvcue_sig <- d$disc_sig
  d$umvcue_sig$beta.exposure <- umvcue(d$disc_sig$beta.exposure, d$rep_sig$beta.exposure, d$disc_sig$se.exposure, d$rep_sig$se.exposure, 5e-8)
  m <- lapply(names(d), function(x) {
    y <- d[[x]]
    if(nrow(y) > 0)
    {
      z <- mr(y, method_list=c("mr_ivw", "mr_wald_ratio")) %>% as_tibble()
      z$what <- x
      return(z)
    } else {
      return(NULL)
    }
  }) %>% bind_rows()
  return(list(bias=bias, mr=m))
}


res <- mclapply(1:1000, function(i) {
  message(i)
  bind_rows(
    gwas_sim(bxy=0.0, out_index = 1:3000, buy=-0.4)$mr %>% mutate(overlap="disc_overlap", bxy=0.0, buy=-0.4),
    gwas_sim(bxy=0.0, out_index = 3001:6000, buy=-0.4)$mr %>% mutate(overlap="rep_overlap", bxy=0.0, buy=-0.4),
    gwas_sim(bxy=0.0, out_index = 6001:9000, buy=-0.4)$mr %>% mutate(overlap="no_overlap", bxy=0.0, buy=-0.4),
    gwas_sim(bxy=0.0, out_index = 1501:4500, buy=-0.4)$mr %>% mutate(overlap="both_overlap", bxy=0.0, buy=-0.4),
    gwas_sim(bxy=0.0, out_index = 1:3000, buy=0)$mr %>% mutate(overlap="disc_overlap", bxy=0.0, buy=0),
    gwas_sim(bxy=0.0, out_index = 3001:6000, buy=0)$mr %>% mutate(overlap="rep_overlap", bxy=0.0, buy=0),
    gwas_sim(bxy=0.0, out_index = 6001:9000, buy=0)$mr %>% mutate(overlap="no_overlap", bxy=0.0, buy=0),
    gwas_sim(bxy=0.0, out_index = 1501:4500, buy=0)$mr %>% mutate(overlap="both_overlap", bxy=0.0, buy=0),
    gwas_sim(bxy=0.0, out_index = 1:3000, buy=0.4)$mr %>% mutate(overlap="disc_overlap", bxy=0.0, buy=0.4),
    gwas_sim(bxy=0.0, out_index = 3001:6000, buy=0.4)$mr %>% mutate(overlap="rep_overlap", bxy=0.0, buy=0.4),
    gwas_sim(bxy=0.0, out_index = 6001:9000, buy=0.4)$mr %>% mutate(overlap="no_overlap", bxy=0.0, buy=0.4),
    gwas_sim(bxy=0.0, out_index = 1501:4500, buy=0.4)$mr %>% mutate(overlap="both_overlap", bxy=0.0, buy=0.4),
    gwas_sim(bxy=0.2, out_index = 1:3000, buy=-0.4)$mr %>% mutate(overlap="disc_overlap", bxy=0.2, buy=-0.4),
    gwas_sim(bxy=0.2, out_index = 3001:6000, buy=-0.4)$mr %>% mutate(overlap="rep_overlap", bxy=0.2, buy=-0.4),
    gwas_sim(bxy=0.2, out_index = 6001:9000, buy=-0.4)$mr %>% mutate(overlap="no_overlap", bxy=0.2, buy=-0.4),
    gwas_sim(bxy=0.2, out_index = 1501:4500, buy=-0.4)$mr %>% mutate(overlap="both_overlap", bxy=0.2, buy=-0.4),
    gwas_sim(bxy=0.2, out_index = 1:3000, buy=0)$mr %>% mutate(overlap="disc_overlap", bxy=0.2, buy=0),
    gwas_sim(bxy=0.2, out_index = 3001:6000, buy=0)$mr %>% mutate(overlap="rep_overlap", bxy=0.2, buy=0),
    gwas_sim(bxy=0.2, out_index = 6001:9000, buy=0)$mr %>% mutate(overlap="no_overlap", bxy=0.2, buy=0),
    gwas_sim(bxy=0.2, out_index = 1501:4500, buy=0)$mr %>% mutate(overlap="both_overlap", bxy=0.2, buy=0),
    gwas_sim(bxy=0.2, out_index = 1:3000, buy=0.4)$mr %>% mutate(overlap="disc_overlap", bxy=0.2, buy=0.4),
    gwas_sim(bxy=0.2, out_index = 3001:6000, buy=0.4)$mr %>% mutate(overlap="rep_overlap", bxy=0.2, buy=0.4),
    gwas_sim(bxy=0.2, out_index = 6001:9000, buy=0.4)$mr %>% mutate(overlap="no_overlap", bxy=0.2, buy=0.4),
    gwas_sim(bxy=0.2, out_index = 1501:4500, buy=0.4)$mr %>% mutate(overlap="both_overlap", bxy=0.2, buy=0.4)
  )
}, mc.cores=16) %>% bind_rows()

save(res, file=file.path(datadir, "replication.rdata"))

