#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")

library(data.table)
require(tidyverse)

args <- commandArgs(T)

datadir <- args[1]
resultsdir <- args[2]
resultsdir <- args[3]

# Read all phenotype names and define each phenotype id
phen_all <- read.table(paste(datadir,"/ukb-b-idlist.txt",sep=""))

# Replace all capital letters with lowercase
phen_all <- phen_all %>% mutate(V1 = tolower(V1))

out_mr <- c()
out_het <- c()

k <- 0

stdevs <- fread(paste(sddir,"ukb-b-sd.csv",sep="/"), header=TRUE)

for (id in phen_all[,1])
{
	# Reading the results for each trait
 	mr_file <- paste0(resultsdir, "/", id, "/MR_All_vs_All.txt")
  het_file <- paste0(resultsdir, "/", id, "/MR_Het_All_vs_All.txt")

  mr_table <- fread(mr_file, header=TRUE)
  het_table <- fread(het_file, header=TRUE)

  mr_table <- mr_table[complete.cases(mr_table[ , "b"]), ]
  het_table <- het_table[complete.cases(het_table[ , "Q"]), ]

  mr_table <- merge(mr_table, stdevs, by.x = "id.exposure", by.y = "id")
  colnames(mr_table)[which(names(mr_table) == "sd")] <- "sd_exp"

  otc=id
  sd_out <- subset(stdevs,stdevs$id == otc)
  sd_list <- rep(sd_out$sd,nrow(mr_table))

  mr_table <- cbind(mr_table,"sd_out"=sd_list)

  mr_table <- cbind(mr_table,"b_tr"=mr_table$b * mr_table$sd_exp / mr_table$sd_out)
  mr_table <- cbind(mr_table,"se_tr"=mr_table$se * mr_table$sd_exp / mr_table$sd_out)

  out_mr <- rbind(out_mr, mr_table)
  out_het <- rbind(out_het, het_table)

  k <- k+1

  print(round(k/length(phen_all[,1])*100))

}

# Remove this later
out_mr <- fread(paste(datadir,"/MR_All_vs_All_combined.txt",sep=""), header=TRUE)

out_mr <- as.data.frame(out_mr)
out_het <- as.data.frame(out_het)

write.table(out_mr, file = paste(datadir,"/MR_All_vs_All_combined.txt",sep=""), append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

write.table(out_het, file = paste(datadir,"/MR_Het_All_vs_All_combined.txt",sep=""), append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


# Plot correlation figures
t_mr <- out_mr
t_mr$pair <- paste(t_mr$id.exposure,t_mr$id.outcome,sep='_')



res_plot <- function(dr,dr_dt,dt,met) {
  t_mr_dat <- t_mr[t_mr$dir==dr_dt&t_mr$data==dt,]
  t_mr_DR <- t_mr[t_mr$dir==dr&t_mr$data=="DR",]
  x <- t_mr_dat[,c("method","b_tr","pair")]
  y <- t_mr_DR[,c("method","b_tr","pair")]
  
  x_met <- subset(x,x$method==met)
  y_met <- subset(y,y$method==met)
  m_met <- merge(x_met, y_met, by.x = "pair", by.y = "pair")
  m_met <- m_met[,c("pair","b_tr.x","b_tr.y")]
  row.names(m_met) <- m_met[,"pair"]
  m_met <- m_met[,-1]
  colnames(m_met) <- c(dt,"DR")

  png(paste(datadir,"/",dr,"_DR_vs_",dr_dt,"_",dt,"_",met,".png",sep=""))
  plot(m_met[,dt], m_met$DR,main="Correlation BETA",xlab=paste(dr_dt,dt,sep="_"),ylab=paste(dr,"DR",sep="_"),xlim=c(min(m_met$D, m_met$DR), max(m_met$D, m_met$DR)), ylim=c(min(m_met$D, m_met$DR), max(m_met$D, m_met$DR)))
  abline(lm(m_met$DR ~ m_met[,dt]))
  abline(coef = c(0,1),col="red")
  dev.off()

  # Plot Histogram
  #png(paste(datadir,"/Hist_",dr,"_DR_vs_",dt,"_",met,".png",sep=""))
  #hist(m_met$D,breaks=200)
  #dev.off()

}



# Results and plots for different tests
res_plot("AB","AB","D","MR Egger")
res_plot("AB","AB","D","Weighted median")
res_plot("AB","AB","D","Inverse variance weighted")
res_plot("AB","AB","D","Simple mode")
res_plot("AB","AB","D","Weighted mode")

res_plot("AB","AB","R","MR Egger")
res_plot("AB","AB","R","Weighted median")
res_plot("AB","AB","R","Inverse variance weighted")
res_plot("AB","AB","R","Simple mode")
res_plot("AB","AB","R","Weighted mode")

res_plot("BA","BA","D","MR Egger")
res_plot("BA","BA","D","Weighted median")
res_plot("BA","BA","D","Inverse variance weighted")
res_plot("BA","BA","D","Simple mode")
res_plot("BA","BA","D","Weighted mode")

res_plot("BA","BA","R","MR Egger")
res_plot("BA","BA","R","Weighted median")
res_plot("BA","BA","R","Inverse variance weighted")
res_plot("BA","BA","R","Simple mode")
res_plot("BA","BA","R","Weighted mode")

res_plot("AB","BA","R","MR Egger")
res_plot("AB","BA","R","Weighted median")
res_plot("AB","BA","R","Inverse variance weighted")
res_plot("AB","BA","R","Simple mode")
res_plot("AB","BA","R","Weighted mode")
