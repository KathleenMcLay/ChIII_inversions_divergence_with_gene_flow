#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=180GB
#SBATCH --job-name=GEA_filter
#SBATCH --time=2:00:00
#SBATCH --partition=general
#SBATCH --account=a_ortiz_barrientos_coe
#SBATCH --output=/home/uqkmcla4/scripts/GEA/GEA_out.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=k.mclay@uq.edu.au

#Filter the VCF to just samples from populations with environmental data 
module load bcftools
module load vcftools
module load htslib
module load intel
module load r

scratch="/scratch/user/uqkmcla4"

bcftools view --threads 12 -O z -S /scratch/user/uqkmcla4/GEA_samples.txt -o ${scratch}/sf8_final_GEA_pops.vcf.gz ${scratch}/sf8_no_variants_noD1_reheader_final_nsng.vcf.gz
bcftools index --threads 12 -t ${scratch}/sf8_final_GEA_pops.vcf.gz

#Filter for MAF 0.03 (should also remove any monomorphic sites because MAF is zero) & max missing 
vcftools --gzvcf ${scratch}/sf8_final_GEA_pops.vcf.gz --maf 0.03 --recode --recode-INFO-all --out ${scratch}/sf8_final_GEA_pops_filtered_MAF_03
vcftools --vcf ${scratch}/sf8_final_GEA_pops_filtered_MAF_03.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out ${scratch}/sf8_final_GEA_pops_filtered_MAF_missing_03

# Prune the dataset to get 10kb putatively neutral SNPs/unlinked SNPs
bcftools +prune -w 195kb --nsites-per-win 1 ${scratch}/sf8_final_GEA_pops_filtered_MAF_missing_03.recode.vcf -O v -o ${scratch}/sf8_final_GEA_pops_filtered_10kb_SNPs_03.vcf 
bcftools stats ${scratch}/sf8_final_GEA_pops_filtered_10kb_SNPs_03.vcf > ${scratch}/03_stats_10kb.txt

# Generate input files from VCF
cat ${scratch}/sf8_final_GEA_pops_filtered_10kb_SNPs_03.vcf | perl /home/uqkmcla4/scripts/reformat/vcf2baypass.pl ${scratch}/GEA_pops.txt ${scratch}/sf8_final_GEA_pops_filtered_10kb_SNPs_03.baypass
cat ${scratch}/sf8_final_GEA_pops_filtered_MAF_missing_03.recode.vcf | perl /home/uqkmcla4/scripts/reformat/vcf2baypass.pl ${scratch}/GEA_pops.txt ${scratch}/sf8_final_GEA_pops_filtered_SNPs_03.baypass

# Run BayPass under the core model mode to generate covariance matrix
/home/uqkmcla4/i_baypass -npop 12 -gfile ${scratch}/sf8_final_GEA_pops_10kb.baypass -outprefix ${scratch}/popcov_ -nthreads 20

# Run BayPass under the standard covariate model
/home/uqkmcla4/i_baypass -npop 12 -gfile ${scratch}/sf8_final_GEA_pops_filtered_SNPs_03.baypass -efile ${scratch}/GEA_soil_data.txt -omegafile ${scratch}/popcov_mat_omega.out -outprefix ${scratch}/all_soilGEA -nthreads 20

# Produce POD samples with 10kb SNP data to determine significance threshold
Rscript /home/uqkmcla4/scripts/ch2/GEA/sim_pop_data.R
/home/uqkmcla4/i_baypass -npop 12 -gfile ${scratch}/G.sf8_final_GEA_pops_10kb.pod -efile ${scratch}/GEA_soil_data.txt -omegafile ${scratch}/popcov_mat_omega.out -outprefix ${scratch}/pod_soilGEA -nthreads 20