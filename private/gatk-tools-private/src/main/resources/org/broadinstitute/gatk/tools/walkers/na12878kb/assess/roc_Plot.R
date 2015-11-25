options(echo=TRUE) 
args <- commandArgs(trailingOnly = TRUE)
print(args)
## args[1] is the report with all the data
## args[2] is the file-name for the ROC-count plots
## args[3] is the file-name for the ROC-ratio plot

library(gsalib)
library(ggplot2)
roc <- gsa.read.gatkreport(args[1])$NA12878Assessment

project <- levels(roc$project)
cplot_fname <- args[2]
rplot_fname <- args[3]

## Following the java code we use labels for index 1 and 2 to keep track of what they represent.
snp_idx   <- 1
indel_idx <- 2

## Variant labels
var_lab <- rep("a", 2)
var_lab[snp_idx]   <- "SNPs"
var_lab[indel_idx] <- "Indels"

## Uncalled row indexes (should be 1 and 2, unless java code has changed).
u_idx<- c(0,0)
u_idx[snp_idx]   <- which(roc$variation == "uncalled_SNPs")
u_idx[indel_idx] <- which(roc$variation == "uncalled_Indels")

## Re-form the uncalled and totalcalled matrices from the java code
uncalled <- matrix(0, nrow=2, ncol=2)
colnames(uncalled) <- c("N_TP", "N_FP")
rownames(uncalled) <- var_lab
uncalled["SNPs",   "N_TP"] <- roc[u_idx[snp_idx],]$N_TP
uncalled["Indels", "N_TP"] <- roc[u_idx[indel_idx],]$N_TP
uncalled["SNPs",   "N_FP"] <- roc[u_idx[snp_idx],]$N_FP
uncalled["Indels", "N_FP"] <- roc[u_idx[indel_idx],]$N_FP

totalcalled <- matrix(0,nrow=2, ncol=2)
colnames(totalcalled) <- c("N_TP", "N_FP")
rownames(totalcalled) <- var_lab
totalcalled["SNPs",   "N_TP"] <- max(c(subset(roc, variation=="SNPs")$N_TP,0))
totalcalled["Indels", "N_TP"] <- max(c(subset(roc, variation=="Indels")$N_TP,0))
totalcalled["SNPs",   "N_FP"] <- max(c(subset(roc, variation=="SNPs")$N_FP,0))
totalcalled["Indels", "N_FP"] <- max(c(subset(roc, variation=="Indels")$N_FP,0))

## Generate the ROC-count plots for SNPs and Indels and save them in the file with name given as arguments.
## Firstly save the appropriate x and y limits in a data frame
TP_lims <- totalcalled[,"N_TP"] + uncalled[,"N_TP"]
FP_lims <- totalcalled[,"N_FP"]
FTP_lim_df <- data.frame(variation=var_lab, TP_lim=TP_lims, FP_lim=FP_lims)

ggp_roc_count <- ggplot(data=roc[-c(u_idx),], aes(x=N_FP, y=N_TP, colour=vqslod)) +
geom_line() +
    geom_hline(data=FTP_lim_df, aes(yintercept=TP_lim), color="gray30", linetype="dashed") +
    geom_text(data=FTP_lim_df, aes(x=FP_lim/2, y=TP_lim + 2*TP_lim/100, label="Total TP in NA12878-DB"), size=3, color="gray30") +
    xlab("Number of False Positives
           Note: line is grey where vqslod is less than 0") +
    ylab("Number of True Positives") +
    facet_wrap(~variation, scales="free") +
    ggtitle("ROC-count Curves") +
    scale_colour_gradient(limits = c(0,max(roc$vqslod)))

ggsave(filename=cplot_fname, plot=ggp_roc_count, type="cairo-png", width=8, height=5, units="in", dpi=200)

## Change the counts to Ratios: TP <- TP/total(TP) etc.
## Firstly for SNps
idx <- which(roc$variation == "SNPs")
roc[idx,]$N_TP <- roc[idx,]$N_TP / totalcalled["SNPs", "N_TP"]
idx <- which(roc$variation == "SNPs")
roc[idx,]$N_FP <- roc[idx,]$N_FP / totalcalled["SNPs", "N_FP"]

## Then for Indels
idx <- which(roc$variation == "Indels")
roc[idx,]$N_TP <- roc[idx,]$N_TP / totalcalled["Indels", "N_TP"]
idx <- which(roc$variation == "Indels")
roc[idx,]$N_FP <- roc[idx,]$N_FP / totalcalled["Indels", "N_FP"]

## Plots the ROC-ratio plot and save it in the file with name given as argument.
ggp_roc_ratio <- ggplot(data=roc[-c(u_idx),], aes(x=N_FP, y=N_TP, linetype=variation)) +
  geom_line() +
  scale_linetype_manual(values=c("dashed", "solid")) +
  xlab("Ratio of False Positives") +
  ylab("Ratio of True Positives") +
  ggtitle("ROC-ratio Curves")
ggsave(filename=rplot_fname, plot=ggp_roc_ratio, type="cairo-png", width=6, height=5, units="in", dpi=200)

