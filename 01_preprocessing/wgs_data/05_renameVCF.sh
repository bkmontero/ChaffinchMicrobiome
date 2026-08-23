#!/bin/bash
#SBATCH -c 4
#SBATCH -n 2 
#SBATCH --mem-per-cpu=2G 
#SBATCH -t 4:00:00

module load samtools

# rename sample names vcf
# rename doesnt work with compressed files
#cd /mnt/netapp1/Store_csddepla/FringillaProject/genotyping_pipeline/output/03-variants_filtered/microbiome_subset
#zcat variants_microbiome_subset_miss0.5_mac3_noZ.vcf.gz | wc -c
cd /mnt/netapp1/Store_csddepla/FringillaProject/genotyping_pipeline/output/03-variants_filtered/microbiome_subset

bcftools view variants_microbiome_subset_miss0.5_mac3_noZ.vcf.gz -O z -o variants_microbiome_subset_miss0.5_mac3_noZ_bgz.vcf.gz

bcftools reheader -s rename_microbiome.tsv -o variants_microbiome_subset_miss0.5_mac3_noZ_renamed.vcf.gz variants_microbiome_subset_miss0.5_mac3_noZ_bgz.vcf.gz

