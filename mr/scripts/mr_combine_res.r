#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")

library(data.table)
require(tidyverse)

args <- commandArgs(T)

datadir <- args[1]
resultsdir <- args[2]
sddir <- args[3]
imgdir <- args[4]

# Read all phenotype names and define each phenotype id
phen_all <- read.table(paste(datadir,"/ukb-b-idlist.txt",sep=""))

# Replace all capital letters with lowercase
phen_all <- phen_all %>% mutate(V1 = tolower(V1))

out_mr <- c()
out_het <- c()

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

  mr_table <- cbind(mr_table,"z"=abs(mr_table$b/mr_table$se))

  mr_table <- mr_table[,-c("V1")]
  
  out_mr <- rbind(out_mr, mr_table)
  out_het <- rbind(out_het, het_table)

}

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
  t_mr_DR <- t_mr[t_mr$dir==dr&t_mr$data=="DR",]
  t_mr_dat <- t_mr[t_mr$dir==dr_dt&t_mr$data==dt,]
  
  x <- t_mr_DR[,c("method","b_tr","z","pair")]
  y <- t_mr_dat[,c("method","b_tr","z","pair")]
  
  x_met <- subset(x,x$method==met)
  y_met <- subset(y,y$method==met)
  m_met <- merge(x_met, y_met, by.x = "pair", by.y = "pair")

  m_met <- m_met[,c("pair","b_tr.x","b_tr.y","z.x","z.y")]
  row.names(m_met) <- m_met[,"pair"]
  m_met <- m_met[,-1]
  colnames(m_met) <- c("DR",dt,"z_DR",paste("z_",dt,sep=""))

  m_met <- m_met[complete.cases(m_met),]

  png(paste(imgdir,"/",dr,"_DR_vs_",dr_dt,"_",dt,"_",met,".png",sep=""))
  plot(m_met$DR,m_met[,dt],main="Correlation BETA",xlab=paste(dr,"DR",sep="_"),ylab=paste(dr_dt,dt,sep="_"), xlim=c(min(m_met[,dt], m_met$DR), max(m_met[,dt], m_met$DR)),ylim=c(min(m_met[,dt], m_met$DR), max(m_met[,dt], m_met$DR)))
  abline(lm(m_met[,dt] ~ m_met$DR))
  abline(coef = c(0,1),col="red")
  dev.off()

  png(paste(imgdir,"/Z_",dr,"_DR_vs_",dr_dt,"_",dt,"_",met,".png",sep=""))
  plot(m_met$z_DR, m_met[,paste("z",dt,sep="_")], main="Correlation Z",xlab=paste("z",dr,"DR",sep="_"),ylab=paste("z",dr_dt,dt,sep="_"),xlim=c(min(m_met[,paste("z",dt,sep="_")], m_met$z_DR), max(m_met[,paste("z",dt,sep="_")], m_met$z_DR)), ylim=c(min(m_met[,paste("z",dt,sep="_")], m_met$z_DR), max(m_met[,paste("z",dt,sep="_")], m_met$z_DR)))
  abline(lm(m_met[,paste("z",dt,sep="_")] ~ m_met$z_DR))
  abline(coef = c(0,1),col="red")
  dev.off()

  # Plot Histogram
  #png(paste(datadir,"/Hist_",dr,"_DR_vs_",dt,"_",met,".png",sep=""))
  #hist(m_met$D,breaks=200)
  #dev.off()

  q <- summary(lm(m_met[,dt] ~ m_met$DR))

  rt <- data.frame("X"=paste(dr,"DR",sep="_"),"Y"=paste(dr_dt,dt,sep="_"),"Method"=met,"Regr_est"=q$coefficients["m_met$DR","Estimate"],"Std_err"=q$coefficients["m_met$DR","Std. Error"])

  return(rt)
}



# Results and plots for different tests
a <- c()
a <- rbind(a,res_plot("AB","AB","D","MR Egger"))
a <- rbind(a,res_plot("AB","AB","D","Weighted median"))
a <- rbind(a,res_plot("AB","AB","D","Inverse variance weighted"))
a <- rbind(a,res_plot("AB","AB","D","Simple mode"))
a <- rbind(a,res_plot("AB","AB","D","Weighted mode"))

a <- rbind(a,res_plot("AB","AB","R","MR Egger"))
a <- rbind(a,res_plot("AB","AB","R","Weighted median"))
a <- rbind(a,res_plot("AB","AB","R","Inverse variance weighted"))
a <- rbind(a,res_plot("AB","AB","R","Simple mode"))
a <- rbind(a,res_plot("AB","AB","R","Weighted mode"))

a <- rbind(a,res_plot("AB","AB","RD","MR Egger"))
a <- rbind(a,res_plot("AB","AB","RD","Weighted median"))
a <- rbind(a,res_plot("AB","AB","RD","Inverse variance weighted"))
a <- rbind(a,res_plot("AB","AB","RD","Simple mode"))
a <- rbind(a,res_plot("AB","AB","RD","Weighted mode"))

a <- rbind(a,res_plot("BA","BA","D","MR Egger"))
a <- rbind(a,res_plot("BA","BA","D","Weighted median"))
a <- rbind(a,res_plot("BA","BA","D","Inverse variance weighted"))
a <- rbind(a,res_plot("BA","BA","D","Simple mode"))
a <- rbind(a,res_plot("BA","BA","D","Weighted mode"))

a <- rbind(a,res_plot("BA","BA","R","MR Egger"))
a <- rbind(a,res_plot("BA","BA","R","Weighted median"))
a <- rbind(a,res_plot("BA","BA","R","Inverse variance weighted"))
a <- rbind(a,res_plot("BA","BA","R","Simple mode"))
a <- rbind(a,res_plot("BA","BA","R","Weighted mode"))

a <- rbind(a,res_plot("BA","BA","RD","MR Egger"))
a <- rbind(a,res_plot("BA","BA","RD","Weighted median"))
a <- rbind(a,res_plot("BA","BA","RD","Inverse variance weighted"))
a <- rbind(a,res_plot("BA","BA","RD","Simple mode"))
a <- rbind(a,res_plot("BA","BA","RD","Weighted mode"))

a <- rbind(a,res_plot("AB","BA","R","MR Egger"))
a <- rbind(a,res_plot("AB","BA","R","Weighted median"))
a <- rbind(a,res_plot("AB","BA","R","Inverse variance weighted"))
a <- rbind(a,res_plot("AB","BA","R","Simple mode"))
a <- rbind(a,res_plot("AB","BA","R","Weighted mode"))

a <- rbind(a,res_plot("AB","BA","RD","MR Egger"))
a <- rbind(a,res_plot("AB","BA","RD","Weighted median"))
a <- rbind(a,res_plot("AB","BA","RD","Inverse variance weighted"))
a <- rbind(a,res_plot("AB","BA","RD","Simple mode"))
a <- rbind(a,res_plot("AB","BA","RD","Weighted mode"))

a <- rbind(a,res_plot("BA","AB","R","MR Egger"))
a <- rbind(a,res_plot("BA","AB","R","Weighted median"))
a <- rbind(a,res_plot("BA","AB","R","Inverse variance weighted"))
a <- rbind(a,res_plot("BA","AB","R","Simple mode"))
a <- rbind(a,res_plot("BA","AB","R","Weighted mode"))

a <- rbind(a,res_plot("BA","AB","RD","MR Egger"))
a <- rbind(a,res_plot("BA","AB","RD","Weighted median"))
a <- rbind(a,res_plot("BA","AB","RD","Inverse variance weighted"))
a <- rbind(a,res_plot("BA","AB","RD","Simple mode"))
a <- rbind(a,res_plot("BA","AB","RD","Weighted mode"))

a <- a[!a$Method=="MR Egger",]

# Save all results
write.table(a, file = paste(datadir,"/Regression_estimates.txt",sep=""), append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
