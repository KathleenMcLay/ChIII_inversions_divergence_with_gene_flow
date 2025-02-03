#!/bin/bash

module load vcftools

DIR="/QRISdata/Q6656/chapter_III/population_divergence"
gzvcf="/QRISdata/Q6656/sf8_no_variants_noD1_reheader_final_nsng.vcf.gz"

# Read population file containing pop name, P1 samples, P2 samples
POP_FILE="/home/uqkmcla4/ChIII_inversions_divergence_with_gene_flow/population_pair_divergence/pops.txt"

# Loop through each population in the first column of the file
while IFS=',' read -r pop P1 P2; do
    # Define output file paths
    FST_OUT="${DIR}/${pop}_fst"
    
    # Create temporary files for sample lists
    echo "$P1" | tr ";" "\n" > "${DIR}/${pop}_P1.txt"
    echo "$P2" | tr ";" "\n" > "${DIR}/${pop}_P2.txt"

    # Calculate FST
    vcftools \
        --gzvcf ${gzvcf} \
        --weir-fst-pop "${DIR}/${pop}_P1.txt" \
        --weir-fst-pop "${DIR}/${pop}_P2.txt" \
        --fst-window-size 10000 --fst-window-step 10000 \
        --out ${FST_OUT}

done < "$POP_FILE"
