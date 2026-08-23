#!/bin/bash
#SBATCH -c 24
#SBATCH -n 2 
#SBATCH --mem-per-cpu=2G 
#SBATCH -t 8:00:00
#SBATCH --mail-type=begin #Envía un correo cuando el trabajo inicia
#SBATCH --mail-type=end #Envía un correo cuando el trabajo finaliza
#SBATCH --mail-user=b.karina.montero@gmail.com #Dirección a la que se envía

module load samtools
module load gcc/system openmpi/4.0.5_ft3 vcflib/1.0.5

CHAFF=/mnt/netapp1/Store_IMIB/FringillaProject/genotyping_pipeline/output//03-variants_filtered/microbiome_subset/variants_microbiome_subset.vcf.gz
OUT=/mnt/netapp1/Store_IMIB/FringillaProject/genotyping_pipeline/output//03-variants_filtered/microbiome_subset/variants_microbiome_subset_miss0.5_mac3_noZ.vcf.gz


# set filters
MAFA=0.05
MAFB=0.1
MISS=0.5
QUAL=30
MIN=3
MAX_DEPTH=16
MAC=3


# perform the filtering with vcftools

vcftools --gzvcf $CHAFF \
		--minQ $QUAL \
		--minDP $MIN \
		--mac $MAC \
        --max-missing $MISS \
        --not-chr OY740764,OY740727 \
		--recode --stdout | gzip -c > \
		$OUT

