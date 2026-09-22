#!/bin/bash
#SBATCH --job-name=redo_vcf
#SBATCH --output=redo_vcf.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH -p general

source /home/trf0019/miniforge3/bin/activate variant_env

REF="/scratch/acicola_project/fastas/sydowiaspp/sydowiafasta/goodbusco/BS06_240_394-4.scaffolds.fa"
BAMS="/scratch/acicola_project/fastas/sydowiaspp/sydowiafastq/Alignmentgood/*.sorted.bam"
OUT="/scratch/acicola_project/fastas/sydowiaspp/sydowiafastq/Alignmentgood/vcf_results/merged_variants_fixed.vcf.gz"

bcftools mpileup -Ou -f "$REF" $BAMS | \
  bcftools call -mv -Oz -o "$OUT"

tabix -p vcf "$OUT"

bcftools query -l "$OUT"
