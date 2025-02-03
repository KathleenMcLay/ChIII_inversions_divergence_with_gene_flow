#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=240GB
#SBATCH --job-name=phy
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --account=a_ortiz_barrientos_coe
#SBATCH --output=/home/uqkmcla4/scripts/phy_out.txt

module load bcftools
module load python
module load vcftools

scratch="/scratch/user/uqkmcla4"

# filter data to just populations with gene flow estimates. 
bcftools view --threads 12 -O z -S ${scratch}/phy/phy_inds.txt -o ${scratch}/phy/sf8_phy_pops.vcf.gz ${scratch}/sf8_no_variants_noD1_reheader_final_nsng.vcf.gz
bcftools index --threads 12 -t ${scratch}/phy/sf8_phy_pops.vcf.gz

# filter for MAF 0.05 and max missing 20%
vcftools --gzvcf ${scratch}/phy/sf8_phy_pops.vcf.gz --maf 0.05 --recode --recode-INFO-all --out ${scratch}/phy/sf8_phy_pops_MAF05
vcftools --vcf ${scratch}/phy/sf8_phy_pops_MAF05.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out ${scratch}/phy/sf8_phy_pops_MAF05_80M_final

# filter to 1 SNP per 10kb 
bcftools +prune -w 10kb --nsites-per-win 1 ${scratch}/phy/sf8_phy_pops_MAF05_80M_final.recode.vcf -O v -o ${scratch}/phy/sf8_phy_pops_MAF05_80M_final_10kb.vcf

# convert to phy file 
/home/uqkmcla4/vcf2phylip.py --input ${scratch}/phy/sf8_phy_pops_MAF05_80M_final_10kb.vcf --output-folder ${scratch}/phy

# RAxML phylogeny 
/home/uqkmcla4/raxmlHPC-PTHREADS-AVX2 -T 46 -f a -d -m GTRGAMMA -p 45687 -N 1000 -e 0.01 -x 56043 -B 0.03 --bootstop-perms 100 -n ML -w ${scratch}/phy -s ${scratch}/phy/sf8_phy_pops_MAF05_80M_final_10kb.min4.phy \
    -o A1003,A1004,A1005,A1008,A1010,A1011
