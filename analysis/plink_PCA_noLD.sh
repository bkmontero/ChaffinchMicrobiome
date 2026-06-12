module load gcc/system
module load plink/2.00a2.3



VCF=/mnt/netapp1/Store_IMIB/FringillaProject/genotyping_pipeline/output/03-variants_filtered/microbiome_subset/variants_microbiome_subset_miss0.8_mac3_noZ.vcf.gz

FINCHES=variants_microbiome_subset_miss0.8_mac3_noZ


# PCA with no LD-pruning
FINCHES=variants_microbiome_subset_miss0.8_mac3_noZ_noLD

# PCA
plink2 --vcf $VCF --double-id --allow-extra-chr \
      --set-missing-var-ids @:# \
      --make-bed --pca --out $FINCHES

