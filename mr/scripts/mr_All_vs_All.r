#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(gwasglue)
library(gwasvcf)
library(ieugwasr)
library(genetics.binaRies)
set_bcftools()

args <- commandArgs(T)

datadir <- args[1]
resultsdir <- args[2]
vcfdir <- args[3]
out <- args[4]


# ======================== Functions ==============================

mr_analysis <- function(exp,d,r,out_gwas) {
  if (!out_gwas=="NA") {
    print(paste("exposure=",exp))
    
    # Read the results of GWAS
    df <- paste(resultsdir,"/",exp,"/",d,".statsfile.txt.gz",sep="")
    dsc <- fread(df,header=TRUE)
    dsc <- as.data.frame(dsc)
    
    rf <- paste(resultsdir,"/",exp,"/",r,".statsfile.txt.gz",sep="")
    rpc <- fread(rf,header=TRUE)
    rpc <- as.data.frame(rpc)
    
    if ("BETA" %in% colnames(dsc) & "BETA" %in% colnames(rpc)) {
      
      dsc1 <- dsc[,c("SNP","ALLELE1","ALLELE0","A1FREQ","BETA", "SE",tail(colnames(dsc),1))]
      colnames(dsc1) <- c("SNP",
                          "effect_allele.exposure",
                          "other_allele.exposure",
                          "eaf.exposure",
                          "beta.exposure",
                          "se.exposure",
                          "pval.exposure")
      dsc2 <- dsc1[order(dsc1$SNP, dsc1$pval.exposure), ]
      dsc2 <- dsc2[ !duplicated(dsc2$SNP), ]  
      dsc3 <- cbind(dsc2,
                    "exposure"=rep("exposure",nrow(dsc2)),
                    "mr_keep.exposure"=rep("TRUE", nrow(dsc2)),
                    "pval_origin.exposure"=rep("reported", nrow(dsc2)),
                    "id.exposure"=rep(exp, nrow(dsc2)),
                    "data_source.exposure"=rep("textfile", nrow(dsc2)))
      
      
      
      
      rpc1 <- rpc[,c("SNP","ALLELE1","ALLELE0","A1FREQ","BETA", "SE",tail(colnames(rpc),1))]
      colnames(rpc1) <- c("SNP",
                          "effect_allele.exposure",
                          "other_allele.exposure",
                          "eaf.exposure",
                          "beta.exposure",
                          "se.exposure",
                          "pval.exposure")
      rpc2 <- rpc1[order(rpc1$SNP, rpc1$pval.exposure), ]
      rpc2 <- rpc2[ !duplicated(rpc2$SNP), ]  
      rpc3 <- cbind(rpc2,
                    "exposure"=rep("exposure",nrow(rpc2)),
                    "mr_keep.exposure"=rep("TRUE", nrow(rpc2)),
                    "pval_origin.exposure"=rep("reported", nrow(rpc2)),
                    "id.exposure"=rep(exp, nrow(rpc2)),
                    "data_source.exposure"=rep("textfile", nrow(rpc2)))
      
      
      disc_gwas <- dsc3
      repl_gwas <- rpc3
      
      
      # Filtering SNPs for their presence in the phenotype and for p-val
      snp_exp <- as.character(subset(dat, id.exposure == exp)$SNP)
      
      disc_gwas_f <- subset(disc_gwas,SNP %in% snp_exp & pval.exposure < 5e-8)
      
      repl_gwas_s <- subset(repl_gwas,SNP %in% snp_exp & pval.exposure < 5e-8)
      
      repl_gwas_f <- subset(repl_gwas,SNP %in% disc_gwas_f$SNP)
      
      out_gwas_f <- subset(out_gwas,SNP %in% disc_gwas_f$SNP)
      
      out_gwas_s <- subset(out_gwas,SNP %in% repl_gwas_s$SNP)
      
      # DR and D scenarios
      if (!nrow(disc_gwas_f)==0) {
        # Harmonise the exposure and outcome data
        dat_dr <- harmonise_data(repl_gwas_f, out_gwas_f, action=1)
        dat_d <- harmonise_data(disc_gwas_f, out_gwas_f, action=1)
        
        # Perform MR on the replication data
        res_dr <- mr(dat_dr)
        if (nrow(res_dr)==0) {
          res_dr <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                               "exposure"="NA","method"="NA","nsnp"=NA,"b"=NA,"se"=NA,"pval"=NA)
        }
        res_dr <- cbind(res_dr,"data"=rep("DR", nrow(res_dr)))
        
        het_dr <- mr_heterogeneity(dat_dr)
        if (nrow(het_dr)==0) {
          het_dr <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                               "exposure"="NA","method"="NA","Q"=NA,"Q_df"=NA,"Q_pval"=NA)
        }
        het_dr <- cbind(het_dr,"data"=rep("DR", nrow(het_dr)))
        
        # Perform MR on the discovery sign data
        res_d <- mr(dat_d)
        if (nrow(res_d)==0) {
          res_d <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                              "exposure"="NA","method"="NA","nsnp"=NA,"b"=NA,"se"=NA,"pval"=NA)
        }
        res_d <- cbind(res_d,"data"=rep("D", nrow(res_d)))
        
        het_d <- mr_heterogeneity(dat_d)
        if (nrow(het_d)==0) {
          het_d <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                              "exposure"="NA","method"="NA","Q"=NA,"Q_df"=NA,"Q_pval"=NA)
        }
        het_d <- cbind(het_d,"data"=rep("D", nrow(het_d)))
        
      } else {
        res_dr <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                             "exposure"="NA","method"="NA","nsnp"=NA,"b"=NA,"se"=NA,"pval"=NA)
        res_dr <- cbind(res_dr,"data"=rep("DR", nrow(res_dr)))
        
        het_dr <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                             "exposure"="NA","method"="NA","Q"=NA,"Q_df"=NA,"Q_pval"=NA)
        het_dr <- cbind(het_dr,"data"=rep("DR", nrow(het_dr)))
        
        res_d <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                            "exposure"="NA","method"="NA","nsnp"=NA,"b"=NA,"se"=NA,"pval"=NA)
        res_d <- cbind(res_d,"data"=rep("D", nrow(res_d)))
        
        het_d <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                            "exposure"="NA","method"="NA","Q"=NA,"Q_df"=NA,"Q_pval"=NA)
        het_d <- cbind(het_d,"data"=rep("D", nrow(het_d)))
      }
      
      
      # R scenario
      if (!nrow(repl_gwas_s)==0) {
        # Harmonise the exposure and outcome data
        dat_r <- harmonise_data(repl_gwas_s, out_gwas_s, action=1)
        
        # Perform MR on the replication sign data
        res_r <- mr(dat_r)
        if (nrow(res_r)==0) {
          res_r <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                              "exposure"="NA","method"="NA","nsnp"=NA,"b"=NA,"se"=NA,"pval"=NA)
        }
        res_r <- cbind(res_r,"data"=rep("R", nrow(res_r)))
        
        het_r <- mr_heterogeneity(dat_r)
        if (nrow(het_r)==0) {
          het_r <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                              "exposure"="NA","method"="NA","Q"=NA,"Q_df"=NA,"Q_pval"=NA)
        }
        het_r <- cbind(het_r,"data"=rep("R", nrow(het_r)))
        
      } else {
        res_r <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                            "exposure"="NA","method"="NA","nsnp"=NA,"b"=NA,"se"=NA,"pval"=NA)
        res_r <- cbind(res_r,"data"=rep("R", nrow(res_r)))
        
        het_r <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                            "exposure"="NA","method"="NA","Q"=NA,"Q_df"=NA,"Q_pval"=NA)
        het_r <- cbind(het_r,"data"=rep("R", nrow(het_r)))
      }
      
    } else {
      res_dr <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                           "exposure"="NA","method"="NA","nsnp"=NA,"b"=NA,"se"=NA,"pval"=NA)
      res_dr <- cbind(res_dr,"data"=rep("DR", nrow(res_dr)))
      
      het_dr <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                           "exposure"="NA","method"="NA","Q"=NA,"Q_df"=NA,"Q_pval"=NA)
      het_dr <- cbind(het_dr,"data"=rep("DR", nrow(het_dr)))
      
      res_d <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                          "exposure"="NA","method"="NA","nsnp"=NA,"b"=NA,"se"=NA,"pval"=NA)
      res_d <- cbind(res_d,"data"=rep("D", nrow(res_d)))
      
      het_d <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                          "exposure"="NA","method"="NA","Q"=NA,"Q_df"=NA,"Q_pval"=NA)
      het_d <- cbind(het_d,"data"=rep("D", nrow(het_d)))
      
      res_r <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                          "exposure"="NA","method"="NA","nsnp"=NA,"b"=NA,"se"=NA,"pval"=NA)
      res_r <- cbind(res_r,"data"=rep("R", nrow(res_r)))
      
      het_r <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                          "exposure"="NA","method"="NA","Q"=NA,"Q_df"=NA,"Q_pval"=NA)
      het_r <- cbind(het_r,"data"=rep("R", nrow(het_r)))
    }
    
  } else {
    res_dr <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                         "exposure"="NA","method"="NA","nsnp"=NA,"b"=NA,"se"=NA,"pval"=NA)
    res_dr <- cbind(res_dr,"data"=rep("DR", nrow(res_dr)))
    
    het_dr <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                         "exposure"="NA","method"="NA","Q"=NA,"Q_df"=NA,"Q_pval"=NA)
    het_dr <- cbind(het_dr,"data"=rep("DR", nrow(het_dr)))
    
    res_d <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                        "exposure"="NA","method"="NA","nsnp"=NA,"b"=NA,"se"=NA,"pval"=NA)
    res_d <- cbind(res_d,"data"=rep("D", nrow(res_d)))
    
    het_d <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                        "exposure"="NA","method"="NA","Q"=NA,"Q_df"=NA,"Q_pval"=NA)
    het_d <- cbind(het_d,"data"=rep("D", nrow(het_d)))
    
    res_r <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                        "exposure"="NA","method"="NA","nsnp"=NA,"b"=NA,"se"=NA,"pval"=NA)
    res_r <- cbind(res_r,"data"=rep("R", nrow(res_r)))
    
    het_r <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                        "exposure"="NA","method"="NA","Q"=NA,"Q_df"=NA,"Q_pval"=NA)
    het_r <- cbind(het_r,"data"=rep("R", nrow(het_r)))
  }
  
  mr_res <- rbind(res_dr,res_d,res_r)
  mr_het <- rbind(het_dr,het_d,het_r)
  
  res_list <- list("mr" = mr_res, "het" = mr_het)
  
  return(res_list)
  
}




outcome <- function(x) {
  of <- paste(resultsdir,"/",out,"/",x,".statsfile.txt.gz",sep="")
  opt <- fread(of,header=TRUE)
  opt <- as.data.frame(opt)
  
  
  if ("BETA" %in% colnames(opt)) {
    opt1 <- opt[,c("SNP","ALLELE1","ALLELE0","A1FREQ","BETA", "SE",tail(colnames(opt),1))]
    colnames(opt1) <- c("SNP",
                        "effect_allele.outcome",
                        "other_allele.outcome",
                        "eaf.outcome",
                        "beta.outcome",
                        "se.outcome",
                        "pval.outcome")
    opt2 <- opt1[order(opt1$SNP, opt1$pval.outcome), ]
    opt2 <- opt2[ !duplicated(opt2$SNP), ]  
    opt3 <- cbind(opt2,
                  "outcome"=rep("outcome",nrow(opt2)),
                  "mr_keep.outcome"=rep("TRUE", nrow(opt2)),
                  "pval_origin.outcome"=rep("reported", nrow(opt2)),
                  "id.outcome"=rep(out, nrow(opt2)),
                  "data_source.outcome"=rep("textfile", nrow(opt2)))
    
    out_gwas <- opt3
  } else {
    out_gwas <- "NA"
  }
  
  return(out_gwas)
}


# ======================== Run ==============================

# Read all phenotype names ands define each phenotype id
phen_all <- fread(paste(datadir,"ukb-b-idlist.txt",sep="/"), header=FALSE)
phen_all <- as.data.frame(phen_all)

# ---------------------- Full MR ----------------------------

# Read instruments for all of the exposures
load(paste(datadir,"MR_prep.RData",sep="/"))
exposure_dat <- mybiglist$exp_dat
chrpos <- mybiglist$chr_pos

# Lookup from one dataset
filename <- paste(vcfdir,"/",out,"/",out,".vcf.gz",sep="")

out1 <- query_gwas(filename, chrompos=chrpos)

out2 <- gwasglue::gwasvcf_to_TwoSampleMR(out1, "outcome")

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, out2, action=1)

# testing for one exposure only: 
# dat <- dat[dat$id.exposure==exp,]

# Perform full MR
res <- mr(dat)
if (nrow(res)==0) {
  res <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                    "exposure"="NA","method"="NA","nsnp"=NA,"b"=NA,"se"=NA,"pval"=NA)
}
res <- res %>% mutate(outcome = tolower(outcome))
res$id.outcome <- res$outcome


het <- mr_heterogeneity(dat)
if (nrow(het)==0) {
  het <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                    "exposure"="NA","method"="NA","Q"=NA,"Q_df"=NA,"Q_pval"=NA)
}
het <- het %>% mutate(outcome = tolower(outcome))
het$id.outcome <- het$outcome


# ---------------------- Replication -------------------------
out_gwas_d <- outcome("discovery")

out_gwas_r <- outcome("replication")

mr_out <- c()
het_out <- c()

for (exp in phen_all[,1])
{ 
  mr_full <- subset(res, id.exposure==exp)
  if (nrow(mr_full)==0) {
    mr_full <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                          "exposure"="NA","method"="NA","nsnp"=NA,"b"=NA,"se"=NA,"pval"=NA)
  }
  mr_full <- cbind(mr_full,"data"=rep("full", nrow(mr_full)))
  mr_full <- cbind(mr_full,"dir"="NA")
  
  het_full <- subset(het, id.exposure==exp)
  if (nrow(het_full)==0) {
    het_full <- data.frame("id.exposure"=exp,"id.outcome"=out,"outcome"="NA",
                           "exposure"="NA","method"="NA","Q"=NA,"Q_df"=NA,"Q_pval"=NA)
  }
  het_full <- cbind(het_full,"data"=rep("full", nrow(het_full)))
  het_full <- cbind(het_full,"dir"="NA")
  
  mr_rep_AB <- mr_analysis(exp,"discovery","replication",out_gwas_r)
  
  mr_AB <- mr_rep_AB$mr
  mr_AB <- cbind(mr_AB,"dir"=rep("AB", nrow(mr_AB)))
  
  het_AB <- mr_rep_AB$het
  het_AB <- cbind(het_AB,"dir"=rep("AB", nrow(het_AB)))
  
  
  
  mr_rep_BA <- mr_analysis(exp,"replication","discovery",out_gwas_d)
  
  mr_BA <- mr_rep_BA$mr
  mr_BA <- cbind(mr_BA,"dir"=rep("BA", nrow(mr_BA)))
  
  het_BA <- mr_rep_BA$het
  het_BA <- cbind(het_BA,"dir"=rep("BA", nrow(het_BA)))
  
  
  mr_out <- rbind(mr_out,mr_full,mr_AB,mr_BA)
  het_out <- rbind(het_out,het_full,het_AB,het_BA)
  
}

# Save all results
write.table(mr_out, file = paste(resultsdir,out,"MR_All_vs_All.txt",sep="/"), append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

write.table(het_out, file = paste(resultsdir,out,"MR_Het_All_vs_All.txt",sep="/"), append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

