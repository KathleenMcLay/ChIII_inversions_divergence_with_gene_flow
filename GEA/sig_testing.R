library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(gridExtra)
library(grid)

############ ANALYSIS OF BAYPASS OUTPUT ############

### PART 1: Process Data ###

# Read inversion data
inversion_data <- read.csv("GEA_inversions_consecutive.csv")
# add a column that combines the first three columns
inversion_data <- inversion_data %>%
  mutate(inversion = paste(scaffold, start, end, sep = "_")) %>%
  mutate(con_start = start / 1000000, con_end = end / 1000000) # convert to Mbp for plotting

# Read in the GEA results file with BF.dB values and scaffold bp position for each SNP
gea_dat <- read.table("positions_all_soilGEA_summary_betai_reg.out", header = TRUE)
gea_dat$position <- gea_dat$position / 1000000
gea_dat$scaffold <- str_replace(gea_dat$scaffold, "scaffold_([0-9]+)", function(x) {
  sprintf("scaffold_%02d", as.integer(str_extract(x, "[0-9]+")))
})
gea_dat <- gea_dat[order(gea_dat$COVARIABLE, gea_dat$scaffold, gea_dat$position),] #reorder the SNPs by scaffold and position

# Adjust the SNP positions so the results can be plotted consecutively for all scaffolds in one plot

snps <- subset(gea_dat, COVARIABLE == 1)[, c("position", "scaffold")] # get a list of SNP positions 
bp <- snps$position[snps$scaffold == "scaffold_01"] # make a new variable with all the positions from scaffold 1
mid <- 0+max(bp)/2 # create a value that is the midpoint of scaffold 1
scaffolds <- unique(snps$scaffold)

for (i in scaffolds[-1]) {
  max.now <- max(bp)
  mid <- c(mid, max.now+max(snps$position[snps$scaffold == i])/2) # mid is for plotting, doesn't get added to gea_dat
  bp <- c(bp, (snps$position[snps$scaffold == i] + max.now))
}
snps$con_pos <- bp #add the new positions to the snps dataframe

# Add consecutive snp position to gea_dat
gea_list <- list()
for (i in unique(gea_dat$COVARIABLE)) {
  gea_dat2 <- gea_dat %>% 
    filter(COVARIABLE == i) %>% 
    inner_join(snps, by = c("scaffold", "position"))  
  
  gea_list[[as.character(i)]] <- gea_dat2
}
gea_dat <- bind_rows(gea_list)

# Flag the inversion status of each SNP
gea_dat <- gea_dat %>%
  rowwise() %>%
  mutate(
    # Find all inversions this SNP falls within
    inversion_name = {
      inv_matches <- which(position >= inversion_data$start & 
                           position <= inversion_data$end & 
                           scaffold == inversion_data$scaffold)
      if (length(inv_matches) > 0) {
        # Combine all matching inversion names into a comma-separated string
        paste(inversion_data$inversion[inv_matches], collapse = ",")
      } else {
        NA_character_  # NA if not in any inversion
      }
    },
    # Set inversion status based on whether we found any matches
    inversion_status = if_else(!is.na(inversion_name), "inversion", "non-inversion")
  ) %>%
  ungroup()

### PART 2: Identify significant GEA SNPs ###

# Create a new column that indicates if the SNP is significant for BF.db > 20
gea_dat$significant_BF.dB <- ifelse(gea_dat$BF.dB. > 20, "TRUE", "FALSE")

# Read in pod data and extract 99th percentile value of Bf.db scores for each covariable
pod_dat <- read.table("pod_soilGEA_summary_betai_reg.out", header=T)
pod_dat_thresh <- pod_dat %>% group_by(COVARIABLE) %>% summarise(pod_threshold = quantile(BF.dB., 0.99))

# Create a new column in gea_dat that indicates if the SNP is significant for 99th percentile of POD
gea_dat <- gea_dat %>%
  left_join(pod_dat_thresh, by = "COVARIABLE") %>%
  mutate(significant_99th_POD = BF.dB. >= pod_threshold)
gea_dat <- gea_dat %>% select(-pod_threshold) #remove the threshold column 

### PART 3: Identify if significant GEA SNPs fall within inversions ###

# Get significant SNP counts per inversion for each covariable 
significant_snp_counts <- inversion_data %>%
  left_join(gea_dat, by = "scaffold", relationship = "many-to-many") %>%
  filter(position >= start & position <= end & significant_BF.dB == TRUE) %>%
  group_by(across(all_of(names(inversion_data))), COVARIABLE) %>%
  summarise(snp_count = n(), .groups = "drop") %>%
  pivot_wider(names_from = COVARIABLE, 
              values_from = snp_count, 
              values_fill = 0,
              names_prefix = "snp_count_cov_")

# Get total SNP counts for each inversion
total_snp_counts <- inversion_data %>%
  left_join(gea_dat, by = "scaffold", relationship = "many-to-many") %>%
  filter(COVARIABLE == "1") %>%  
  filter(position >= start & position <= end) %>%
  group_by(across(all_of(names(inversion_data)))) %>%
  summarise(total_snps = n(), .groups = "drop")

# Combine 
snp_counts <- significant_snp_counts %>%
  left_join(total_snp_counts, by = names(inversion_data))

# For each inversions Calculate the proportion of SNPs that are significant each covariable
snp_counts <- snp_counts %>%
  mutate(across(starts_with("snp_count_cov_"), ~ .x / total_snps, .names = "prop_{col}"))


# Calculate the proportion of significant SNPs in each inversion of the total significant SNPs for that covariable using the gea_summary data and snp_counts
snp_counts <- snp_counts %>%
  mutate(
    prop_cov_1_sig = snp_count_cov_1 / gea_summary$num_SNPs_BF.dB[gea_summary$COVARIABLE == 1],
    prop_cov_2_sig = snp_count_cov_2 / gea_summary$num_SNPs_BF.dB[gea_summary$COVARIABLE == 2],
    prop_cov_3_sig = snp_count_cov_3 / gea_summary$num_SNPs_BF.dB[gea_summary$COVARIABLE == 3],
)

### PART 4: Create a summary of the data ###





# Create a summary of the data 
gea_summary <- gea_dat %>%
  left_join(pod_dat_thresh, by = "COVARIABLE") %>%
  left_join(gea_dat %>% select(scaffold, position, COVARIABLE, inversion_status), 
            by = c("scaffold", "position", "COVARIABLE")) %>%
  group_by(COVARIABLE) %>%
  summarise(
    # Significant GEA SNPs counts
    tot_sig_BF.dB = sum(significant_BF.dB == TRUE, na.rm = TRUE),
    tot_sig_POD.99 = sum(significant_99th_POD == TRUE, na.rm = TRUE),
    sig_POD.99 = first(pod_threshold),
    
    # Significant GEA SNPs in inversions
    count_BF.dB_inv = sum(significant_BF.dB == TRUE & inversion_status == "inversion", na.rm = TRUE),
    count_POD.99_inv = sum(significant_99th_POD == TRUE & inversion_status == "inversion", na.rm = TRUE),
    prop_BF.dB_inv = count_BF.dB_inv / tot_sig_BF.dB,
    prop_POD.99_inv = count_POD.99_inv / tot_sig_POD.99,
    
    # Total SNPs in inversions  
    tot_SNPs = n(),
    tot_SNPs_invs = sum(inversion_status == "inversion", na.rm = TRUE),
    prop_tot_SNPs_invs = tot_SNPs_invs / n(),
    .groups = "drop"
  )

### PART 5: Run permutation testing: random sampling of significant # SNPs within inversions ###

# Initialize an empty tibble to store final permutation results
prop_inv_perm <- tibble()

# Set number of permutation iterations
num_iterations <- 1000

# Initialize lists to store results
results_list <- list()  # For summary statistics
counter <- 1  # Counter to keep track of list indices
all_SNPs_out <- list()  # For detailed permutation results

# Loop through each environmental covariable in the GEA results
for (u in unique(gea_summary$COVARIABLE)) {
  
  # Count how many significant SNPs exist for this covariable 
  total_SNPs <- gea_summary$tot_sig_BF.dB[gea_summary$COVARIABLE == u]
  print(paste("Total SNPs for covariable", u, ":", total_SNPs))
  
  # Create a dataframe to store results for each permutation iteration
  SNPs_out <- data.frame(
    iteration = integer(num_iterations),        # Permutation iteration number
    within_inversion = integer(num_iterations), # Count of SNPs within inversions
    perc = numeric(num_iterations),             # Percentage of SNPs within inversions
    COVARIABLE = character(num_iterations),     # Environmental covariable being analyzed
    stringsAsFactors = FALSE
  )
  
  # Run permutation iterations
  for (i in 1:num_iterations) {
    print(paste("permutation no ", i))
    
    # Initialize counter for SNPs found within inversions
    
    snp_count <- 0
    
    # Filter data for the current covariable
    gea_dat_fil <- gea_dat %>% filter(COVARIABLE == u)  
    
    # For each iteration, sample the same number of SNPs as were found significant in the actual analysis
    for (j in 1:total_SNPs) { 
      # Randomly sample a SNP from the entire dataset
      sim_SNP <- sample(1:nrow(gea_dat_fil), 1)
      
      # Get position and scaffold for the randomly sampled SNP
      sim_SNP_position <- gea_dat$position[sim_SNP] 
      sim_SNP_scaffold <- gea_dat$scaffold[sim_SNP]
      
      # Check if the SNP falls within any inversions, first filter to scaffold
      matching_inversions <- inversion_data[inversion_data$scaffold == sim_SNP_scaffold, ]
      
      # If there is an inversion on the scaffold of the current sim_SNP_position, check each inversion to see if the SNP falls within it
      if (nrow(matching_inversions) > 0) {
        for (k in 1:nrow(matching_inversions)) {
          if (sim_SNP_position >= matching_inversions$start[k] && 
              sim_SNP_position <= matching_inversions$end[k]) {
            # If the SNP is within the inversion, increment count and stop checking other inversions
            snp_count <- snp_count + 1
            break  # Break out of the loop once we find one matching inversion
          }
        }
      }
    }
    
    # Store results for this iteration
    SNPs_out[i, ] <- list(i, snp_count, snp_count / total_SNPs, u)
  }
  
  # Store detailed results for each covariable 
  all_SNPs_out[[counter]] <- SNPs_out
  
  # Calculate summary statistics and store in results list
  results_list[[counter]] <- tibble(
    COVARIABLE = u,
    mean_prop_inv_perm = mean(as.numeric(SNPs_out$perc), na.rm = TRUE),  # Average proportion across permutations
    max_prop_inv_perm = max(SNPs_out$perc),                              # Maximum proportion observed
    min_prop_inv_perm = min(SNPs_out$perc)                               # Minimum proportion observed
  )
  
  # Increment counter for next covariable
  counter <- counter + 1
}

# Combine all summary results into a single tibble
prop_inv_perm <- bind_rows(results_list)

# add prop_inv_perm to gea_summary
gea_summary <- gea_summary %>%
  left_join(prop_inv_perm, by = "COVARIABLE") %>%
  mutate(
    mean_prop_inv_perm = ifelse(is.na(mean_prop_inv_perm), 0, mean_prop_inv_perm),
    max_prop_inv_perm = ifelse(is.na(max_prop_inv_perm), 0, max_prop_inv_perm),
    min_prop_inv_perm = ifelse(is.na(min_prop_inv_perm), 0, min_prop_inv_perm)
  )

# Combine all detailed permutation results into a single dataframe
final_SNPs_out <- bind_rows(all_SNPs_out)

### PART 6: Write results to files ###

# Write detailed permutation results to file
write.table(final_SNPs_out, "permutation_SNPs_out.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Write summary results to file 
write.table(gea_summary, "GEA_summary.tsv", sep = "\t", row.names = FALSE)
# Save gea_dat with significant SNPs + inversion status 
write.table(gea_dat, "positions_all_soilGEA_summary_betai_reg_SIG_INV.out", sep = "\t", row.names = FALSE, quote = FALSE)

# Save SNP counts to file
write.table(snp_counts, "per_inversion_sig_SNPS.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

### PART 7: Plotting ###

# Manhatten Plot with GEA data and inversions
for (i in unique(gea_dat$COVARIABLE)) {
  
  plot_data <- gea_dat %>% 
    filter(COVARIABLE == i) %>% 
    select(scaffold, con_pos, BF.dB.) 
  
  p <- ggplot() + 
    theme_classic() +
    geom_point(data = plot_data, mapping = aes(x = con_pos, y = BF.dB., col = scaffold), 
               pch = 20, size = 1, na.rm = TRUE, show.legend = FALSE) +
    geom_hline(aes(yintercept = 20), col = "black", alpha = 0.5, linetype = "dotted") +
    geom_rect(data = inversion_data, inherit.aes = FALSE,
              aes(xmin = con_start, xmax = con_end,
                  ymin = min(21), ymax = max(25)),
              fill = "blue") +
    scale_colour_manual(values = rep(c("lightgrey", "darkgrey"), 7)) +
    scale_y_continuous(expand = c(0, 1), limits = c(-5, 55)) +
    scale_x_continuous(
      breaks = mid,
      labels = substr(scaffolds, nchar(scaffolds) - 1, nchar(scaffolds)),
      expand = c(0.02, 0.02)
    ) +
    theme(
      title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",
      axis.line.x = element_blank()
    )
  
  # Save the plot
  ggsave(
    filename = paste0("GEA_plot_", i, "with_inversions.png"),
    plot = p,
    width = 6,
    height = 1.5,
    dpi = 300,
    units = "in"
  )
}

### PART 8: Allele frequencies ###

# Create a dataframe of just signficant SNPs within inversions 
gea_dat_sig_inv <- gea_dat %>%
  filter(significant_BF.dB == "TRUE" & inversion_status == "inversion")
gea_dat_sig_inv$position <- as.numeric(gea_dat_sig_inv$position)

# Read in the allele frequency data 
allele_data <- read.table("positions_sf8_final_GEA_pops_filtered_SNPs_03.baypass", header = FALSE)
allele_data <- cbind(snps[, c("scaffold", "position", "con_pos")], allele_data) #bind snps to allele data
allele_data$position <- as.numeric(allele_data$position)
colnames(allele_data) <- c("scaffold", "position", "con_pos", "D00_R", "D00_A", "D01_R", "D01_A", "D03_R", "D03_A", "D04_R", "D04_A", "D12_R", "D12_A", "H00_R", "H00_A", "H01_R", "H01_A", "H02_R", "H02_A", "H05_R", "H05_A", "H14_R", "H14_A")

# For each covariable, filter allele_data to just significant SNPs within inversions

# Create a dataframe to hold the output
significant_alleles <- tibble()

for (cov in unique(gea_dat_sig_inv$COVARIABLE)) {
  inv_subset <- gea_dat_sig_inv %>% filter(COVARIABLE == cov)
  
  # Filter allele data for significant SNPs
  allele_subset <- allele_data %>%
    semi_join(inv_subset, by = c("scaffold", "position"))

  #add covariable column
  allele_subset <- allele_subset %>%
    mutate(COVARIABLE = cov)
  
  # bind the results from all three covariables into a new dataframe
  significant_alleles <- bind_rows(significant_alleles, allele_subset)
  
}

pop_ids <- unique(gsub("_[RA]$", "", names(significant_alleles)[c(-1, -2, -3, -24)]))

# For each population, calculate frequency of the alt alelle 
for (pop in pop_ids) {
  ref_col <- paste0(pop, "_R")
  alt_col <- paste0(pop, "_A")
  freq_col <- paste0(pop, "_freq")
  
  significant_alleles[[freq_col]] <- ifelse(
    (significant_alleles[[ref_col]] + significant_alleles[[alt_col]]) == 0,
    NA,
    significant_alleles[[alt_col]] / (significant_alleles[[ref_col]] + significant_alleles[[alt_col]])
  )
}

# Keep only scaffold, position, and freq columns
freq_data <- significant_alleles %>%
  select(scaffold, position, COVARIABLE, ends_with("_freq"))

# Bind the freq data with the dataframe of just signficant SNPs within inversions 
allele_freqs <- gea_dat_sig_inv %>%
  left_join(freq_data, by = c("scaffold", "position", "COVARIABLE"))

# create a dataset with a unique row for each SNP + inversion combination
allele_freqs_2 <- allele_freqs %>%
  # Split the 'inversion_name' column into multiple rows if there are ";" separators
  separate_rows(inversion_name, sep = ",") %>%
  # Remove any leading/trailing spaces after splitting
  mutate(inversion_name = str_trim(inversion_name))

#for each covariable, and unique value in inversion_name, take a weighted mean of the frequencies using the BF.dB value as the weight
weighted_freqs <- allele_freqs_2 %>%
  group_by(COVARIABLE, inversion_name) %>%
  summarise(across(ends_with("_freq"), ~ weighted.mean(.x, BF.dB., na.rm = TRUE), .names = "weighted_{col}")) %>%
  ungroup()

#round the frequencies to 2 decimal places
weighted_freqs <- weighted_freqs %>%
  mutate(across(starts_with("weighted_"), ~ round(.x, 2)))

#write the weighted frequencies to a file
write.table(weighted_freqs, "GEA_sig_inv_weighted_mean_allele_frequencies.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Reshape allele frequencies to long format
allele_freqs_long <- weighted_freqs %>%
  pivot_longer(
    cols = starts_with("weighted_"),
    names_to = "population",
    names_prefix = "weighted_",
    names_transform = list(population = ~ gsub("_freq$", "", .)),
    values_to = "allele_frequency"
  )

# Read in the environmental data 
env_data <- read.csv("/Users/kathleenmclay/Library/CloudStorage/GoogleDrive-mclay.kathleen@gmail.com/My Drive/PhD/ChIII_gene_flow/2_soil_pca/PCA_soil_data_names.csv")
env_data <- env_data[1:3,]
colnames(env_data)[1] <- "covariable"  
env_data$covariable <- gsub("PC", "", env_data$covariable) 
env_data$covariable <- as.integer(env_data$covariable)
env_data <- env_data %>%
  mutate(across(-covariable, ~ round(.x, 2)))
env_data <- env_data %>%
  mutate(across(-covariable, as.numeric))

# Reshape environmental data to long format
env_long <- env_data %>%
  pivot_longer(
    cols = -covariable,
    names_to = "population",
    values_to = "env_value"
  )

# Join datasets on covariable and population
plot_data <- allele_freqs_long %>%
  rename(covariable = COVARIABLE) %>%
  inner_join(env_long, by = c("covariable", "population"))

# Split the data by covariable
plot_data_by_cov <- plot_data %>%
  group_by(covariable) %>%
  group_split()

# Loop through each covariable-specific dataset
for (cov_df in plot_data_by_cov) {
  cov <- unique(cov_df$covariable)
  
  # Generate list of plots for this covariable
  plot_list <- cov_df %>%
    group_by(inversion_name) %>%
    group_split() %>%
    lapply(function(df) {
      inv <- unique(df$inversion_name)
      
      ggplot(df, aes(x = env_value, y = allele_frequency)) +
        geom_point(size = 1.5) +
        ggrepel::geom_text_repel(aes(label = population), size = 2.5, max.overlaps = 10) +
        geom_smooth(method = "lm", se = FALSE, color = "blue", size = 0.6) +
        labs(
          title = paste("Inversion:", inv, "| Covariable:", cov),
          x = "Environmental Value",
          y = "Allele Frequency"
        ) +
        theme_minimal(base_size = 8) +
        theme(
          plot.title = element_text(size = 9),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          plot.margin = margin(5, 5, 5, 5)
        )
    })
  
  # Save the plot list into a PDF named by covariable
  pdf_filename <- paste0("allele_freq_vs_environment_", cov, ".pdf")
  pdf(pdf_filename, width = 11, height = 14)
  gridExtra::marrangeGrob(grobs = plot_list, ncol = 3, nrow = 5, top = NULL) %>% print()
  dev.off()
}