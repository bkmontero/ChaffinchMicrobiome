#!/bin/bash
#SBATCH -c 4
#SBATCH -n 2
#SBATCH --mem-per-cpu 2G 
#SBATCH -t 06:00:00

module load samtools

## compute_Fst

# set environm variables
VCF=/mnt/netapp1/Store_IMIB/FringillaProject/genotyping_pipeline/output/03-variants_filtered/microbiome_subset/variants_microbiome_subset_miss0.5_mac3_noZ_renamed.vcf.gz
OUT=/mnt/netapp1/Store_IMIB/FringillaProject/09_WGS_Workflow/MS_microbiome/microbiome_subset_roh.txt

bcftools roh -G30 --AF-dflt 0.4 $VCF -o $OUT


# Get ROH regions (lines that start with "RG")
FILTER=/mnt/netapp1/Store_IMIB/FringillaProject/09_WGS_Workflow/MS_microbiome/microbiome_subset_roh_rg.txt

grep "^RG" $OUT > $FILTER

head $OUTPUT
wc -l $OUTPUT 
#1974161


REF=/mnt/netapp1/Store_IMIB/FringillaProject/reference/bFriCoe1_complete.fa

#Genome lenght (autosomes)
# Sum only OY7407* scaffolds, excluding Z chromosome (OY740727)
awk '$1 ~ /^OY7407/ && $1 != "OY740727" {sum += $2} END {print "Autosomal genome length: " sum " bp"}' ${REF}.fai
#Autosomal genome length: 1005935441 bp
