#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=180GB
#SBATCH --job-name=GEA_filter
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --account=a_ortiz_barrientos_coe
#SBATCH --output=/home/uqkmcla4/scripts/GEA_ALL.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=k.mclay@uq.edu.au

module load bcftools
module load vcftools
module load htslib
module load intel
module load r

dir="/QRISdata/Q6656/chapter_III/GEA/all/final"
scratch="/scratch/user/uqkmcla4"

#Filter the VCF to just samples from populations with environmental data 
bcftools view --threads 12 -O z -S ${dir}/GEA_samples.txt -o ${dir}/sf8_final_GEA_pops.vcf.gz ${scratch}/sf8_no_variants_noD1_reheader_final_nsng.vcf.gz
bcftools index --threads 12 -t ${dir}/sf8_final_GEA_pops.vcf.gz

#Filter for MAF 0.03 (should also remove any monomorphic sites because MAF is zero) & max missing 
vcftools --gzvcf ${dir}/sf8_final_GEA_pops.vcf.gz --maf 0.03 --recode --recode-INFO-all --out ${dir}/sf8_final_GEA_pops_filtered_MAF_03
vcftools --vcf ${dir}/sf8_final_GEA_pops_filtered_MAF_03.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out ${dir}/sf8_final_GEA_pops_filtered_MAF_missing_03

# Prune the dataset to get 10kb putatively neutral SNPs/unlinked SNPs
bcftools +prune -w 195kb --nsites-per-win 1 ${dir}/sf8_final_GEA_pops_filtered_MAF_missing_03.recode.vcf -O v -o ${dir}/sf8_final_GEA_pops_filtered_10kb_SNPs_03.vcf 
bcftools stats ${dir}/sf8_final_GEA_pops_filtered_10kb_SNPs_03.vcf > ${dir}/03_stats_10kb.txt

# Generate input files from VCF
cat ${dir}/sf8_final_GEA_pops_filtered_10kb_SNPs_03.vcf | perl /home/uqkmcla4/scripts/reformat/vcf2baypass.pl ${dir}/GEA_pops.txt ${dir}/sf8_final_GEA_pops_filtered_10kb_SNPs_03.baypass
cat ${dir}/sf8_final_GEA_pops_filtered_MAF_missing_03.recode.vcf | perl /home/uqkmcla4/scripts/reformat/vcf2baypass.pl ${dir}/GEA_pops.txt ${dir}/sf8_final_GEA_pops_filtered_SNPs_03.baypass

#Run BayPass under the core model mode to generate covariance matrix
/home/uqkmcla4/i_baypass -npop 10 -gfile ${dir}/sf8_final_GEA_pops_filtered_10kb_SNPs_03.baypass -outprefix ${dir}/popcov -nthreads 20

#Run BayPass under the standard covariate model
/home/uqkmcla4/i_baypass -npop 10 -gfile ${dir}/sf8_final_GEA_pops_filtered_SNPs_03.baypass -efile ${dir}/PCA_soil_data.txt -omegafile ${dir}/popcov_mat_omega.out -outprefix ${dir}/all_soilGEA -nthreads 20

# Run BayPass under the standard covariate model - INVERSIONS as biallelic loci 
/home/uqkmcla4/i_baypass -npop 10 -gfile ${dir}/inversion_genotypes.baypass -efile ${dir}/PCA_soil_data.txt -omegafile ${dir}/popcov_mat_omega.out -outprefix ${dir}/inversions_soilGEA -nthreads 20


# Produce POD samples with 10kb SNP data to determine significance threshold
Rscript /home/uqkmcla4/ChIII_inversions_divergence_with_gene_flow/GEA/sim_pop_data_all.R
/home/uqkmcla4/i_baypass -npop 10 -gfile ${dir}/G.sf8_final_GEA_pops_10kb.pod -efile ${dir}/PCA_soil_data.txt -omegafile ${dir}/popcov_mat_omega.out -outprefix ${dir}/pod_soilGEA -nthreads 20

### Add scaffold and position columns to the BayPass output
cd ${dir}
# Extract scaffold and position columns from the VCF file (excluding the header)
awk '!/^#/ {print $1, $2}' OFS="\t" "sf8_final_GEA_pops_filtered_MAF_missing_03.recode.vcf" > "all_chrom.txt"
Duplicate the content three times and add a header
cat "all_chrom.txt" "all_chrom.txt" "all_chrom.txt" > temp.txt
echo -e "scaffold\tposition" | cat - temp.txt > "all_chromall.txt"
rm temp.txt 

# Merge with GEA output
paste "all_chrom.txt" "all_soilGEA_summary_betai_reg.out" > "positions_all_soilGEA_summary_betai_reg.out"

