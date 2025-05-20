library(data.table)
library(mvtnorm)
library(geigen)

source("~/baypass_public-v2.41/utils/baypass_utils.R")

setwd("/QRISdata/Q6656/chapter_III/GEA/all")

# Read in the covariance matrix
omega <- as.matrix(read.table("/QRISdata/Q6656/chapter_III/GEA/all/popcov_mat_omega.out"))
# Get estimates (post. mean) of both the a_pi and b_pi parameters of the Pi Beta distribution
pi_beta_coef <- read.table("/QRISdata/Q6656/chapter_III/GEA/all/popcov_summary_beta_params.out",h=T)$Mean
# dataset 
geno_data <- geno2YN("/QRISdata/Q6656/chapter_III/GEA/all/sf8_final_GEA_pops_filtered_10kb_SNPs_03.baypass")
# Create the POD
pod_data <- simulate.baypass(omega.mat=omega, nsnp=10000, sample.size=geno_data$NN, 
                            beta.pi=pi_beta_coef, pi.maf=0, suffix="sf8_final_GEA_pops_10kb.pod")

