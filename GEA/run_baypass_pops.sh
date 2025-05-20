#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=180GB
#SBATCH --job-name=GEA
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --account=a_ortiz_barrientos_coe
#SBATCH --output=/home/uqkmcla4/scripts/GEA_INVERSIONS.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=k.mclay@uq.edu.au

# Load the required modules 
module load bcftools
module load vcftools
module load htslib
module load intel
module load r

#set -x  # Print each command as it executes
set -e  # Exit on any error

#scratch="/scratch/user/uqkmcla4"
directory="/QRISdata/Q6656/chapter_III/GEA"
#home="/home/uqkmcla4/scripts/reformat"

# for file in ${directory}/sample_lists/*; do
#     echo "current file is ${file}"
#     filename=$(basename ${file} _samples.txt)
    #Filter the VCF to just the current population
    # bcftools view --threads 12 -O z \
    #     -S ${file} \
    #     -o ${directory}/${filename}.vcf.gz \
    #     ${scratch}/sf8_no_variants_noD1_reheader_final_nsng.vcf.gz
    # bcftools index --threads 12 -t ${directory}/${filename}.vcf.gz

    # # Filter for MAF 0.03 (should also remove any monomorphic sites because MAF is zero) & max missing 20%
    # vcftools --gzvcf ${directory}/${filename}.vcf.gz \
    #     --maf 0.03 \
    #     --recode --recode-INFO-all --out ${directory}/${filename}_filtered_MAF_03
    # vcftools --vcf ${directory}/${filename}_filtered_MAF_03.recode.vcf \
    #     --max-missing 0.8 \
    #     --recode --recode-INFO-all --out ${directory}/${filename}_filtered_MAF_03_80

    # # Prune the dataset to get 10kb putatively neutral SNPs/unlinked SNPs
    # bcftools +prune -w 195kb --nsites-per-win 1 -O v \
    #     -o ${directory}/${filename}_filtered_SNPs_03_80_10kb.vcf \
    #     ${directory}/${filename}_filtered_MAF_03_80.recode.vcf
    # bcftools stats ${directory}/${filename}_filtered_SNPs_03_80_10kb.vcf > ${directory}/${filename}_stats_10kb.txt

    # Generate input files from VCF
#     cat ${directory}/${filename}_filtered_SNPs_03_80_10kb.vcf | perl ${home}/vcf2baypass.pl ${directory}/${filename}_pops.tsv ${directory}/${filename}_filtered_SNPs_03_80_10kb.baypass
#     cat ${directory}/${filename}_filtered_MAF_03_80.recode.vcf | perl ${home}/vcf2baypass.pl ${directory}/${filename}_pops.tsv ${directory}/${filename}_filtered_SNPs_03_80.baypass

#     # Run BayPass under the core model mode to generate covariance matrix
#     /home/uqkmcla4/i_baypass -npop 2 -gfile ${directory}/${filename}_filtered_SNPs_03_80_10kb.baypass -outprefix ${directory}/${filename}_popcov -nthreads 20

#     # Run BayPass under the standard covariate model
#     /home/uqkmcla4/i_baypass -npop 2 -gfile ${directory}/${filename}_filtered_SNPs_03_80.baypass -efile ${directory}/${filename}_PCA_soil_data.txt -omegafile ${directory}/${filename}_popcov_mat_omega.out -outprefix ${directory}/${filename}_soilGEA -nthreads 20

#     # # Produce POD samples with 10kb SNP data to determine significance threshold
#     Rscript /home/uqkmcla4/ChIII_inversions_divergence_with_gene_flow/GEA/sim_pop_data.R ${filename}
#     /home/uqkmcla4/i_baypass -npop 2 -gfile ${directory}/G.${filename}_filtered_SNPs_03_80_10kb.pod -efile ${directory}/${filename}_PCA_soil_data.txt -omegafile ${directory}/${filename}_popcov_mat_omega.out -outprefix ${directory}/${filename}_pod_soilGEA -nthreads 20

#     cd ${directory}
#     # Extract scaffold and position columns from the VCF file (excluding the header)
#     awk '!/^#/ {print $1, $2}' OFS="\t" "${filename}_filtered_MAF_03_80.recode.vcf" > "${filename}_chrom.txt"
#     # Duplicate the content three times and add a header
#     cat "${filename}_chrom.txt" "${filename}_chrom.txt" "${filename}_chrom.txt" > temp.txt
#     echo -e "scaffold\tposition" | cat - temp.txt > "${filename}_chromall.txt"
#     rm temp.txt 

#     # Merge with GEA output
#     paste "${filename}_chromall.txt" "${filename}_soilGEA_summary_betai_reg.out" > "${filename}_positions_all_soilGEA_summary_betai_reg.out"

# done

# Run baypass on inversions as biallelic SNPs

inv_dir=("/QRISdata/Q6656/chapter_III/GEA/inversions")

array=(
    "D00H00"
    "D01H01"
    "D03H02"
    "D04H05"
    "D12H14"
)

for i in "${array[@]}"; do
    echo "current file is ${i}"
    # Run BayPass under the standard covariate model for inversions 
    /home/uqkmcla4/i_baypass -npop 2 -gfile ${inv_dir}/${i}_inversions.baypass -efile ${directory}/${i}_PCA_soil_data.txt -omegafile ${directory}/${i}_popcov_mat_omega.out -outprefix ${inv_dir}/${i}_soilGEA_inversions -nthreads 20
done

