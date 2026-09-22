#!/usr/bin/env bash
#SBATCH --partition=gpu2,gpu4            # Specify the GPU partition
#SBATCH --gres=gpu:tesla:2               # Request 2 Tesla GPUs
#SBATCH --nodes=1                        # Request 1 node
#SBATCH --ntasks-per-node=20             # 20 tasks per node
#SBATCH --mem=50G                        # Request 50 GB of memory
#SBATCH --time=60:00:00                  # Set the time limit to 10 hours
#SBATCH --job-name=var_assessment      # Job name
#SBATCH --output=filtersnp.out
#SBATCH --error=filtersp.err

# Load necessary modules (if using HPC) or activate environment
source /home/trf0019/miniforge3/bin/activate variant_env  # Uncomment if needed

VCF="/home/trf0019/Working_script/filtered_variants.vcf.gz"
OUT="/scratch/acicola_project/fastas/sydowiaspp/vcf_results/filtered_variants"

# STEP 1 — Keep only SNPs, remove indels, keep bi-allelic
bcftools view \
  -m2 -M2 \
  -v snps \
  -Oz -o ${OUT}_step1.vcf.gz \
  $VCF

tabix -p vcf ${OUT}_step1.vcf.gz

# STEP 2 — Filter low quality and low depth
bcftools filter \
  -e 'QUAL<30 || DP<5' \
  -Oz -o ${OUT}_step2.vcf.gz \
  ${OUT}_step1.vcf.gz

tabix -p vcf ${OUT}_step2.vcf.gz

# STEP 3 — Remove SNPs with missing data (-e filters out missing)
bcftools view \
  -i 'F_MISSING<0.2' \
  -Oz -o ${OUT}.vcf.gz \
  ${OUT}_step2.vcf.gz

tabix -p vcf ${OUT}.vcf.gzVCF="merged_variants.vcf.gz"
OUT="filtered_variants"

# STEP 1 — Keep only SNPs, remove indels, keep bi-allelic
bcftools view \
  -m2 -M2 \
  -v snps \
  -Oz -o ${OUT}_step1.vcf.gz \
  $VCF

tabix -p vcf ${OUT}_step1.vcf.gz

# STEP 2 — Filter low quality and low depth
bcftools filter \
  -e 'QUAL<30 || DP<5' \
  -Oz -o ${OUT}_step2.vcf.gz \
  ${OUT}_step1.vcf.gz

tabix -p vcf ${OUT}_step2.vcf.gz

# STEP 3 — Remove SNPs with missing data (-e filters out missing)
bcftools view \
  -i 'F_MISSING<0.2' \
  -Oz -o ${OUT}.vcf.gz \
  ${OUT}_step2.vcf.gz

tabix -p vcf ${OUT}.vcf.gz