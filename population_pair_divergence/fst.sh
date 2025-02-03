#!/bin/bash

module load vcftools

DIR="/QRISdata/Q6656/chapter_III/population_divergence"
gzvcf="/QRISdata/Q6656/sf8_no_variants_noD1_reheader_final_nsng.vcf.gz"

# Read population file containing P1 samples in column 1 and P2 samples in column 2
POP_FILE="/home/uqkmcla4/ChIII_inversions_divergence_with_gene_flow/population_pair_divergence/pops.txt"

# Loop through each line of the population file
while IFS=$'\t' read -r P1 P2; do
    # Define output file paths
    FST_OUT="${DIR}/${P1}_${P2}_fst"
    
    # Create temporary files for sample lists by filtering sample names starting with P1 and P2
    zgrep -m 1 "^#CHROM" ${gzvcf} | tr '\t' '\n' | grep "^${P1}" > "${DIR}/${P1}_P1.txt"
    zgrep -m 1 "^#CHROM" ${gzvcf} | tr '\t' '\n' | grep "^${P2}" > "${DIR}/${P2}_P2.txt"

    # Calculate FST
    vcftools \
        --gzvcf ${gzvcf} \
        --weir-fst-pop "${DIR}/${P1}_P1.txt" \
        --weir-fst-pop "${DIR}/${P2}_P2.txt" \
        --fst-window-size 10000 --fst-window-step 10000 \
        --out ${FST_OUT}

done < "$POP_FILE"

