library(data.table)
library(mvtnorm)
library(geigen)

source("~/baypass_public-v2.41/utils/baypass_utils.R")

setwd("/QRISdata/Q6656/chapter_III/GEA")

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]  # First argument passed

# Read in the covariance matrix
omega <- as.matrix(read.table(paste(filename, "_popcov_mat_omega.out", sep="")))
# Get estimates (post. mean) of both the a_pi and b_pi parameters of the Pi Beta distribution
pi_beta_coef <- read.table(paste(filename, "_popcov_summary_beta_params.out", sep=""), header=T)$Mean
# dataset 
geno_data <- geno2YN(paste(filename, "_filtered_SNPs_03_80_10kb.baypass", sep=""))
# Create the POD
suffix <- paste(filename, "_filtered_SNPs_03_80_10kb.pod", sep="")
pod_data <- simulate.baypass(omega.mat=omega, nsnp=10000, sample.size=geno_data$NN, 
                            beta.pi=pi_beta_coef, pi.maf=0, suffix=suffix)