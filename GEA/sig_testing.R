library(tidyverse)
library(gridExtra)
library(grid)
select=dplyr::select

# Read in datasets
gea_dat <- read.table("positions_all_soilGEA_summary_betai_reg.out", h=T) # SNP GEA merged with position and scaffold columns of VCF
gea_dat <- gea_dat %>%
  filter(as.numeric(sub("scaffold_", "", scaffold)) < 16) #subset to just first 15 scaffolds 
gea_dat$position <- gea_dat$position/(10^6) #Mb scale

pod_dat <- read.table("pod_soilGEA_summary_betai_reg.out", h=T) # POD GEA for significance threshold
pcs <- data.frame(covariable = c("Soil PC1", "Soil PC2", "Soil PC3")) #names to be used as title for each plot
inv_dat <- read.csv("inversion_list.csv", h=T) # Inversions 
inv_regions <- read.csv("inversion_list_no_overlap.csv", h=T) # inversions regions with overlap removed

#adjust the SNP positions so the results can be plotted consecutively for all scaffolds in one plot
snps <- subset(gea_dat, COVARIABLE == 1)[, c("position", "scaffold")] #get a list of SNP positions 
bp <- snps$position[snps$scaffold == "scaffold_1"] #make a new variable with all the positions from scaffold 1
mid <- 0+max(bp)/2 #create a value that is the midpoint of scaffold 1
scaffolds <- unique(gea_dat$scaffold)
for (i in scaffolds[-1]) {
  print(i)
  max.now <- max(bp)
  mid <- c(mid, max.now+max(snps$position[snps$scaffold == i])/2)
  bp <- c(bp, (snps$position[snps$scaffold == i] + max.now))
  inv_dat$con_start[inv_dat$scaffold == i] = inv_dat$start[inv_dat$scaffold == i] + max.now
  inv_dat$con_end[inv_dat$scaffold == i] = inv_dat$end[inv_dat$scaffold == i] + max.now
}
snps$con_pos <- bp

#plot GEA results for each co-variable + chromosome
pdf("BayPass_plot.pdf", width=12, height=2)
for (i in unique(gea_dat$COVARIABLE)) {
  for (j in unique(gea_dat$scaffold)) {
    plot_dat <- gea_dat %>% 
      filter(COVARIABLE == i, scaffold == j) %>% select(scaffold, position, BF.dB.) #select the current co-variable and chromosome from the GEA data
    my.q <- pod_dat %>%
      filter(COVARIABLE == i) %>% summarize(q99=quantile(BF.dB.,probs=0.99)) %>% pull() # get the significance threshold
    p <- ggplot() + theme_classic() +
      geom_point(data=my.plot, mapping=aes(x=position, y=BF.dB.), col="darkgrey", pch=20, size=2, na.rm=T, show.legend=FALSE) +
      geom_hline(aes(yintercept=my.q), col="#1A242F", alpha=1, linetype="dotted") + # significance level
      geom_hline(aes(yintercept=20), col="#1A242F", alpha=1) + 
      geom_rect(data = subset(inv_dat, scaffold == j), inherit.aes = F,
                aes(xmin = start, xmax = end,
                    ymin = min(0), ymax = max(my.plot$BF.dB.)), fill = "#C1C4AD", alpha = 0.3) +
      xlab("Scaffold") + ylab(expression(BF[is]*" "*(dB))) +
      scale_y_continuous(expand=c(0,0.6), limits=c(0,max(my.plot$BF.dB.))) +
      scale_x_continuous(breaks = seq(0, max(my.plot$position), by = 50)) +
      theme(axis.line.x=element_blank()) +
      ggtitle(paste(pcs$covariable[i], " - ", j))
    print(p)
  }
}
dev.off()

print(p)

#extract significant SNPs for each covariable 
for (i in unique(gea_dat$COVARIABLE)) {
  gea_sig <- gea_dat %>% 
    filter(BF.dB. > 20, COVARIABLE == i) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) 
  assign(paste0("gea_sig_", as.character(i)), gea_sig)
}

# Identify significant SNPs that are within inversions

filtered_gea_sig3 <- data.frame()

for (i in unique(gea_sig_3$scaffold)) {
  scaffold_data <- gea_sig_3 %>% filter(scaffold == i) #filter to current scaffold
  inv_regions_scaffold <- inv_regions %>% filter(scaffold == i) #filter inversions to current scaffold
  filtered_snps <- scaffold_data %>%
    rowwise() %>%
    filter(any(position >= inv_regions_scaffold$start & position <= inv_regions_scaffold$end)) %>%
    ungroup() #check if the SNP is within any of the inversions 
  filtered_gea_sig3 <- bind_rows(filtered_gea_sig3, filtered_snps) #if it is, add it to the new dataset 
}

#total SNPs within inversions 
cov1_total_inv_SNPs <- data.frame()

for (i in unique(gea_dat_cov1$scaffold)) {
  scaffold_data <- gea_dat_cov1 %>% filter(scaffold == i) #filter to current scaffold
  inv_regions_scaffold <- inv_regions %>% filter(scaffold == i) #filter inversions to current scaffold
  inv_snps <- scaffold_data %>%
    rowwise() %>%
    filter(any(position >= inv_regions_scaffold$start & position <= inv_regions_scaffold$end)) %>%
    ungroup() #check if the SNP is within any of the inversions 
  cov1_total_inv_SNPs <- bind_rows(cov1_total_inv_SNPs, inv_snps) #if it is, add it to the new dataset 
}

# Permutation testing - randomise inversions 
inv_reg <- inv_regions[, 5] / 10^6  # Assuming this extracts region lengths
num_iterations <- 1000
perms_out <- data.frame(iteration = integer(num_iterations), sig_BF = integer(num_iterations), total_snps = integer(num_iterations))
gea_dat_cov1 <- gea_dat %>% filter(COVARIABLE == 1)

for (i in 1:num_iterations) {
  
  # Initialize a counter for SNPs with BF > 20
  snp_count <- 0
  total_snps <- 0
  sampled_regions <- logical(length(inv_reg))  
  
  for (j in 1:length(inv_reg)) {
    if (sampled_regions[j]) { # Check if the region has already been sampled
      next  
    }
    sampled_regions[j] <- TRUE # Mark this region as sampled
    region_length <- inv_reg[j] # Get the region length (in terms of position)
    print(paste("Region", j, "Length:", region_length))
    
    # Initialize a flag to track whether a valid end_idx is found
    valid_end_idx_found <- FALSE
    
    # Resample until a valid end_idx is found
    while (!valid_end_idx_found) {
      start_idx <- sample(1:(nrow(gea_dat_cov1) - 1), 1)  
      start_position <- gea_dat_cov1$position[start_idx]
      start_scaffold <- gea_dat_cov1$scaffold[start_idx] 
      end_position <- start_position + region_length
      end_idx <- which(gea_dat_cov1$position >= end_position & gea_dat_cov1$scaffold == start_scaffold)[1]
      
      # Ensure end_idx is after the start_idx
      if (!is.na(end_idx) && end_idx > start_idx) {
        valid_end_idx_found <- TRUE
        region_data <- gea_dat_cov1[start_idx:end_idx, ]
        snp_count <- snp_count + sum(region_data$BF.dB. > 20)
        print(paste("Region", j, "SNPs with BF > 20:", sum(region_data$BF.dB. > 20)))
        total_snps <- total_snps + nrow(region_data)
        print(paste("Region", j, "total SNPs:", nrow(region_data)))
      } else {
        print(paste("Region", j, "could not find a valid end_idx. Resampling..."))
      }
    }
  }
  perms_out[i, ] <- c(i, snp_count, total_snps)
}
# create perc column 
perms_out$perc_snps <- perms_out$sig_BF/perms_out$total_snps

# Permutation testing - randomise significant SNPs  
num_iterations <- 1000
significant_SNPs <- 3516
SNPs_out <- data.frame(iteration = integer(num_iterations), within_inversion = integer(num_iterations))
gea_dat_cov1 <- gea_dat %>% filter(COVARIABLE == 1)

for (i in 1:num_iterations) {
  snp_count <- 0
  for (j in 1:significant_SNPs) { 
    sim_SNP <- sample(1:(nrow(gea_dat_cov1) - 1), 1)
    sim_SNP_position <- gea_dat_cov1$position[sim_SNP] 
    sim_SNP_scaffold <- gea_dat_cov1$scaffold[sim_SNP]
    matching_inversions <- inv_regions[inv_regions$scaffold == sim_SNP_scaffold, ] 
    if (nrow(matching_inversions) > 0) {
      for (k in 1:nrow(matching_inversions)) {
        if (sim_SNP_position >= matching_inversions$start[k] && 
            sim_SNP_position <= matching_inversions$end[k]) {
          snp_count <- snp_count + 1
          break 
        }
      }
    }
  }
  SNPs_out[i, ] <- c(i, snp_count)
}
# create perc column 
SNPs_out$total_snps <- 454311
SNPs_out$perc_snps <- SNPs_out$within_inversion/454311