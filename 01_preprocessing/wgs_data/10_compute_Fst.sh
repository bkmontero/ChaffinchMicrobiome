#!/bin/bash
#SBATCH -c 4
#SBATCH -n 2
#SBATCH --mem-per-cpu 1G 
#SBATCH -t 08:00:00

module load samtools

## compute_Fst

# set environm variables
OUT=/mnt/netapp1/Store_csddepla/FringillaProject/Microbiome/fst
VCF=/mnt/netapp1/Store_csddepla/FringillaProject/genotyping_pipeline/output/03-variants_filtered/microbiome_subset/variants_microbiome_subset_miss0.5_mac3_noZ_renamed.vcf.gz

cd $OUT

# extract sample names per island
bcftools query -l $VCF | grep "F_teydea" > F_teydea
bcftools query -l $VCF | grep "F_polatzeki" > F_polatzeki
bcftools query -l $VCF | grep "Mainland" > Mainland
bcftools query -l $VCF | grep "Sao_Miguel" > Sao_Miguel
bcftools query -l $VCF | grep "Flores" > Flores
bcftools query -l $VCF | grep "Madeira" > Madeira
bcftools query -l $VCF | grep "Hierro" > Hierro
bcftools query -l $VCF | grep "bakeri" > bakeri
bcftools query -l $VCF | grep "Tenerife_F_c_canariensis" > Tenerife_canariensis
bcftools query -l $VCF | grep "La_Gomera_F_c_canariensis" > Gomera_canariensis



vcftools --gzvcf ${VCF} \
--weir-fst-pop Tenerife_canariensis \
--weir-fst-pop Gomera_canariensis \
--out ./Tenerife_canariensis_Gomera_canariensis
#After filtering, kept 21 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.070298
#Weir and Cockerham weighted Fst estimate: 0.13598
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop bakeri \
--weir-fst-pop Tenerife_canariensis \
--out ./bakeri_Tenerife_canariensis
#After filtering, kept 20 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.19697
#Weir and Cockerham weighted Fst estimate: 0.41849
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop bakeri \
--weir-fst-pop Gomera_canariensis \
--out ./bakeri_Gomera_canariensis
#After filtering, kept 21 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.16441
#Weir and Cockerham weighted Fst estimate: 0.36778
#After filtering, kept 29691294 out of a possible 29691294 Sites



vcftools --gzvcf ${VCF} \
--weir-fst-pop F_teydea \
--weir-fst-pop F_polatzeki \
--out ./teydea_polatzeki
#After filtering, kept 20 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.55538
#Weir and Cockerham weighted Fst estimate: 0.79231
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop F_teydea \
--weir-fst-pop Mainland \
--out ./teydea_Mainland
#After filtering, kept 11 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.66679
#Weir and Cockerham weighted Fst estimate: 0.83363
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop F_teydea \
--weir-fst-pop Sao_Miguel \
--out ./teydea_Sao_Miguel
#After filtering, kept 19 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.38442
#Weir and Cockerham weighted Fst estimate: 0.67072
#After filtering, kept 29691294 out of a possible 29691294 Sites


vcftools --gzvcf ${VCF} \
--weir-fst-pop F_teydea \
--weir-fst-pop Flores \
--out ./teydea_Flores
#After filtering, kept 19 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.39613
#Weir and Cockerham weighted Fst estimate: 0.6724
#After filtering, kept 29691294 out of a possible 29691294 Sites


vcftools --gzvcf ${VCF} \
--weir-fst-pop F_teydea \
--weir-fst-pop Madeira \
--out ./teydea_Madeira
#After filtering, kept 20 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.65552
#Weir and Cockerham weighted Fst estimate: 0.8589
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop F_teydea \
--weir-fst-pop Hierro \
--out ./teydea_Hierro
#After filtering, kept 21 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.66474
#Weir and Cockerham weighted Fst estimate: 0.87128
#After filtering, kept 29691294 out of a possible 29691294 Sites


vcftools --gzvcf ${VCF} \
--weir-fst-pop F_teydea \
--weir-fst-pop Palma \
--out ./teydea_Palma
#After filtering, kept 20 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.73821
#Weir and Cockerham weighted Fst estimate: 0.90375
#After filtering, kept 29691294 out of a possible 29691294 Sites


vcftools --gzvcf ${VCF} \
--weir-fst-pop F_teydea \
--weir-fst-pop bakeri \
--out ./teydea_bakeri
#After filtering, kept 19 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.71003
#Weir and Cockerham weighted Fst estimate: 0.89838
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop F_teydea \
--weir-fst-pop Tenerife_canariensis \
--out ./teydea_Tenerife_canariensis
#After filtering, kept 19 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.73935
#Weir and Cockerham weighted Fst estimate: 0.90758
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop F_teydea \
--weir-fst-pop Gomera_canariensis \
--out ./teydea_Gomera_canariensis
#After filtering, kept 20 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.65593
#Weir and Cockerham weighted Fst estimate: 0.87631
#After filtering, kept 29691294 out of a possible 29691294 Sites



vcftools --gzvcf ${VCF} \
--weir-fst-pop F_polatzeki \
--weir-fst-pop Mainland \
--out ./polatzeki_Mainland
#After filtering, kept 13 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.69333
#Weir and Cockerham weighted Fst estimate: 0.85461
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop F_polatzeki \
--weir-fst-pop Sao_Miguel \
--out ./polatzeki_Sao_Miguel
#After filtering, kept 21 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.40189
#Weir and Cockerham weighted Fst estimate: 0.69549
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop F_polatzeki \
--weir-fst-pop Flores \
--out ./polatzeki_Flores
#After filtering, kept 21 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.4134
#Weir and Cockerham weighted Fst estimate: 0.69677
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop F_polatzeki \
--weir-fst-pop Madeira \
--out ./polatzeki_Madeira
#After filtering, kept 22 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.66286
#Weir and Cockerham weighted Fst estimate: 0.87
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop F_polatzeki \
--weir-fst-pop Hierro \
--out ./polatzeki_Hierro
#After filtering, kept 23 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.67096
#Weir and Cockerham weighted Fst estimate: 0.88101
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop F_polatzeki \
--weir-fst-pop Palma \
--out ./polatzeki_Palma
#After filtering, kept 22 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.7418
#Weir and Cockerham weighted Fst estimate: 0.91099
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop F_polatzeki \
--weir-fst-pop bakeri \
--out ./polatzeki_bakeri
#After filtering, kept 21 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.7153
#Weir and Cockerham weighted Fst estimate: 0.90657
#After filtering, kept 29691294 out of a possible 29691294 Sites

#
vcftools --gzvcf ${VCF} \
--weir-fst-pop F_polatzeki \
--weir-fst-pop Tenerife_canariensis \
--out ./polatzeki_Tenerife_canariensis
#After filtering, kept 21 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.74319
#Weir and Cockerham weighted Fst estimate: 0.91477
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop F_polatzeki \
--weir-fst-pop Gomera_canariensis \
--out ./polatzeki_Gomera_canariensis
#After filtering, kept 22 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.66216
#Weir and Cockerham weighted Fst estimate: 0.88581
#After filtering, kept 29691294 out of a possible 29691294 Sites


vcftools --gzvcf ${VCF} \
--weir-fst-pop Mainland \
--weir-fst-pop Sao_Miguel \
--out ./Mainland_Sao_Miguel
#After filtering, kept 12 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.12745
#Weir and Cockerham weighted Fst estimate: 0.22618
#After filtering, kept 29691294 out of a possible 29691294 Sites


vcftools --gzvcf ${VCF} \
--weir-fst-pop Mainland \
--weir-fst-pop Flores \
--out ./Mainland_Flores
#After filtering, kept 12 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.13945
#Weir and Cockerham weighted Fst estimate: 0.23208
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop Mainland \
--weir-fst-pop Madeira \
--out ./Mainland_Madeira
#After filtering, kept 13 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.46302
#Weir and Cockerham weighted Fst estimate: 0.60899
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop Mainland \
--weir-fst-pop Hierro \
--out ./Mainland_Hierro
#After filtering, kept 14 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.48941
#Weir and Cockerham weighted Fst estimate: 0.64453
#After filtering, kept 29691294 out of a possible 29691294 Sites


vcftools --gzvcf ${VCF} \
--weir-fst-pop Mainland \
--weir-fst-pop Palma \
--out ./Mainland_Palma
#After filtering, kept 13 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.54472
#Weir and Cockerham weighted Fst estimate: 0.70493
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop Mainland \
--weir-fst-pop bakeri \
--out ./Mainland_bakeri
#After filtering, kept 12 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.50461
#Weir and Cockerham weighted Fst estimate: 0.67391
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop Mainland \
--weir-fst-pop Tenerife_canariensis \
--out ./Mainland_Tenerife_canariensis
#After filtering, kept 12 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.53298
#Weir and Cockerham weighted Fst estimate: 0.70152
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop Mainland \
--weir-fst-pop Gomera_canariensis \
--out ./Mainland_Gomera_canariensis
#After filtering, kept 13 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.48109
#Weir and Cockerham weighted Fst estimate: 0.65119
#After filtering, kept 29691294 out of a possible 29691294 Sites


vcftools --gzvcf ${VCF} \
--weir-fst-pop Sao_Miguel \
--weir-fst-pop Flores \
--out ./Sao_Miguel_Flores
#After filtering, kept 20 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.01967
#Weir and Cockerham weighted Fst estimate: 0.028991
#After filtering, kept 29691294 out of a possible 29691294 Sites


vcftools --gzvcf ${VCF} \
--weir-fst-pop Sao_Miguel \
--weir-fst-pop Madeira \
--out ./Sao_Miguel_Madeira
#After filtering, kept 21 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.24845
#Weir and Cockerham weighted Fst estimate: 0.44996
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop Sao_Miguel \
--weir-fst-pop Hierro \
--out ./Sao_Miguel_Hierro
#After filtering, kept 22 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.27099
#Weir and Cockerham weighted Fst estimate: 0.4959
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop Sao_Miguel \
--weir-fst-pop Palma \
--out ./Sao_Miguel_Palma
#After filtering, kept 21 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.28477
#Weir and Cockerham weighted Fst estimate: 0.52183
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop Sao_Miguel \
--weir-fst-pop bakeri \
--out ./Sao_Miguel_bakeri
#After filtering, kept 20 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.26661
#Weir and Cockerham weighted Fst estimate: 0.49689
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop Sao_Miguel \
--weir-fst-pop Tenerife_canariensis \
--out ./Sao_Miguel_Tenerife_canariensis
#After filtering, kept 20 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.27776
#Weir and Cockerham weighted Fst estimate: 0.51458
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop Sao_Miguel \
--weir-fst-pop Gomera_canariensis \
--out ./Sao_Miguel_Gomera_canariensis
#After filtering, kept 21 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.26558
#Weir and Cockerham weighted Fst estimate: 0.49358
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop Flores \
--weir-fst-pop Madeira \
--out ./Flores_Madeira
#After filtering, kept 21 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.25959
#Weir and Cockerham weighted Fst estimate: 0.45847
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop Flores \
--weir-fst-pop Hierro \
--out ./Flores_Hierro
#After filtering, kept 22 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.28184
#Weir and Cockerham weighted Fst estimate: 0.50281
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop Flores \
--weir-fst-pop Palma \
--out ./Flores_Palma
#After filtering, kept 21 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.29647
#Weir and Cockerham weighted Fst estimate: 0.52834
#After filtering, kept 29691294 out of a possible 29691294 Sites



vcftools --gzvcf ${VCF} \
--weir-fst-pop Flores \
--weir-fst-pop bakeri \
--out ./Flores_bakeri
#After filtering, kept 20 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.27809
#Weir and Cockerham weighted Fst estimate: 0.50378
#After filtering, kept 29691294 out of a possible 29691294 Sites
#

vcftools --gzvcf ${VCF} \
--weir-fst-pop Flores \
--weir-fst-pop Tenerife_canariensis \
--out ./Flores_Tenerife_canariensis
#After filtering, kept 20 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.28963
#Weir and Cockerham weighted Fst estimate: 0.52127
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop Flores \
--weir-fst-pop Gomera_canariensis \
--out ./Flores_Gomera_canariensis
#After filtering, kept 21 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.27634
#Weir and Cockerham weighted Fst estimate: 0.50058
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop Madeira \
--weir-fst-pop Hierro \
--out ./Madeira_Hierro
#After filtering, kept 23 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.27847
#Weir and Cockerham weighted Fst estimate: 0.52588
#After filtering, kept 29691294 out of a possible 29691294 Sites



vcftools --gzvcf ${VCF} \
--weir-fst-pop Madeira \
--weir-fst-pop Palma \
--out ./Madeira_Palma
#After filtering, kept 22 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.31209
#Weir and Cockerham weighted Fst estimate: 0.55464
#After filtering, kept 29691294 out of a possible 29691294 Sites


vcftools --gzvcf ${VCF} \
--weir-fst-pop Madeira \
--weir-fst-pop bakeri \
--out ./Madeira_bakeri
#After filtering, kept 21 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.29444
#Weir and Cockerham weighted Fst estimate: 0.55432
#After filtering, kept 29691294 out of a possible 29691294 Sites


vcftools --gzvcf ${VCF} \
--weir-fst-pop Madeira \
--weir-fst-pop Tenerife_canariensis \
--out ./Madeira_Tenerife_canariensis
#After filtering, kept 21 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.31531
#Weir and Cockerham weighted Fst estimate: 0.56796
#After filtering, kept 29691294 out of a possible 29691294 Sites


vcftools --gzvcf ${VCF} \
--weir-fst-pop Madeira \
--weir-fst-pop Gomera_canariensis \
--out ./Madeira_Gomera_canariensis
#After filtering, kept 22 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.27093
#Weir and Cockerham weighted Fst estimate: 0.51916
#After filtering, kept 29691294 out of a possible 29691294 Sites


vcftools --gzvcf ${VCF} \
--weir-fst-pop Hierro \
--weir-fst-pop Palma \
--out ./Hierro_Palma
#After filtering, kept 23 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.10628
#Weir and Cockerham weighted Fst estimate: 0.19481
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop Hierro \
--weir-fst-pop bakeri \
--out ./Hierro_bakeri
#After filtering, kept 22 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.20485
#Weir and Cockerham weighted Fst estimate: 0.43583
#After filtering, kept 29691294 out of a possible 29691294 Sites

vcftools --gzvcf ${VCF} \
--weir-fst-pop Hierro \
--weir-fst-pop Tenerife_canariensis \
--out ./Hierro_Tenerife_canariensis
#After filtering, kept 22 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.17888
#Weir and Cockerham weighted Fst estimate: 0.36645
#After filtering, kept 29691294 out of a possible 29691294 Sites


vcftools --gzvcf ${VCF} \
--weir-fst-pop Hierro \
--weir-fst-pop Gomera_canariensis \
--out ./Hierro_Gomera_canariensis
#After filtering, kept 23 out of 107 Individuals
#Outputting Weir and Cockerham Fst estimates.
#Weir and Cockerham mean Fst estimate: 0.1396
#Weir and Cockerham weighted Fst estimate: 0.29129
#After filtering, kept 29691294 out of a possible 29691294 Sites