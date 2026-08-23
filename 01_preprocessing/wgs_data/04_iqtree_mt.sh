#!/bin/bash
#SBATCH -c 8
#SBATCH -n 1 
#SBATCH --mem-per-cpu=1G 
#SBATCH -t 01:00:00


module load cesga/2020 mafft/7.525-with-extensions
module load cesga/2020 iq-tree/2.3.1

# mafft align

FASTA=/mnt/netapp1/Store_csddepla/FringillaProject/nextflow_v2/consensus/mt_Fringilla.fasta
srun mafft --thread 16 $MT/mt_Fringilla.fasta > $MT/mt_Fringilla.fasta_aln.fasta

OUT=/mnt/netapp1/Store_IMIB/FringillaProject/09_WGS_Workflow/MS_microbiome
srun iqtree2 -s $OUT/mt_Fringilla_aln.fasta -m GTR -nt 8 -bb 1000 -seed 80808 -redo