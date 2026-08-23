#!/bin/bash
#SBATCH -c 1
#SBATCH -n 1 
#SBATCH --mem-per-cpu=4G 
#SBATCH -t 99:00:00




### SETTINGS (1 / 2): User input ----------------- ###

    # PROVIDE PATH TO INPUT CSV
    # Note: One row per sample, with unheadered columns "ID, LANE, R1_FASTQ_PATH, R2_FASTQ_PATH"
    sample_csv=/mnt/netapp1/Store_IMIB/FringillaProject/Backup/nextflow_v2/microbiome_subset.csv 
    
    # PROVIDE PATH TO THIS REPOSITORY
    repository_path=/mnt/netapp1/Store_IMIB/FringillaProject/genotyping_pipeline 

    # NAME A NEXTFLOW CONFIGURATION PROFILE
    nextflow_profile="cesga_finisterrae" 

    # SELECT STEPS TO RUN
    # Note: Must be run in order, but filter_variants and phase_variants can be repeated with different settings
    trim_align_reads=yes
    call_variants=yes
    filter_variants=yes
    phase_variants=yes

    # SET OPTIONS FOR TRIM AND ALIGN STEP
    # Note: The read_target is applied separately to R1 and R2 (if downsample_reads=yes)
    # Note: The default flag, 0x400, is for optical and PCR duplicates
    deduplicate_reads=yes
    downsample_reads=yes
    read_target=40000000
    exclude_flags=0x400

    # SET OPTIONS FOR CALL VARIANTS STEP
    concatenate_unfiltered_vcfs=no

    # SET OPTIONS FOR FILTER VARIANTS STEP
    # Note: To refilter output from call_variants, change filtering_label and re-run the filter_variants step
    filtering_label="modern_samples"
    filtering_min_alleles=2
    filtering_max_alleles=2
    filtering_max_missing=0.8
    filtering_min_meanDP=5
    filtering_max_meanDP=16 
    filtering_minDP=5
    filtering_maxDP=16 
    filtering_minQ=30
    filtering_keep="/mnt/netapp1/Store_IMIB/FringillaProject/genotyping_pipeline/modern_samples.tsv"
    
    #filtering_max_meanDP=16 #30 default
    #filtering_maxDP=16 #30 default
    # SET OPTIONS FOR PHASE VARIANTS STEP
    phasing_window_size=10000000

    # PROVIDE DETAILS OF REFERENCE GENOME
    ref_genome=/mnt/netapp1/Store_IMIB/FringillaProject/reference/bFriCoe1_complete.fa
    ref_index=/mnt/netapp1/Store_IMIB/FringillaProject/reference/bFriCoe1_complete.fa.fai # NB: Assuming the index is here, too?
    ref_recombination_map_dir=path_here # NB: This setting is only needed if you want phasing. If so, do you have recombination maps?
    ref_scaffold_name="CAUPSF" # NB: This must be set based on your reference genome
    ref_ploidy_file=./pipeline/assets/default.ploidy # NB: For the time being everything is called as diploid - bugfix incoming within the next few months.

### --------------- End user input --------------- ###

### SETTINGS (2 / 2): Set up environment --------- ###

module load nextflow/25.04.6
#module load cesga/2025 apptainer/1.4.2

### ------------ End set up environment ---------- ###

# Main script begins
set -o errexit
set -o nounset

# Function to handle missing output directories
mkmissingdir() {
    if [ ! -e $1 ]; then
        mkdir -p $1
    fi
}

# Function to check previous step is done
chkprevious() {
    if [ ! -e $2 ]; then
        echo "Error. ${1} expected output from previous step to exist in directory ${2}."
        exit
    fi
}

cd $repository_path
output_dir=${repository_path}/output
trim_align_output_dir=${output_dir}/01-aligned_reads
call_variants_output_dir=${output_dir}/02-variants_unfiltered
filter_variants_output_dir=${output_dir}/03-variants_filtered/${filtering_label}
phase_variants_output_dir=${output_dir}/03-variants_filtered/${filtering_label}/phased
mkmissingdir $output_dir

if [ $trim_align_reads = "yes" ]; then
    mkmissingdir $trim_align_output_dir
    nextflow -log ./.nextflow/nextflow.log \
        run ./pipeline/nextflow/trim_align_reads.nf \
        -with-report $trim_align_output_dir/workflow_report.html \
        -profile $nextflow_profile \
        -resume \
        --samples $sample_csv \
        --deduplicate $deduplicate_reads \
        --downsample $downsample_reads \
        --read_target $read_target \
        --exclude_flags $exclude_flags \
        --ref_genome $ref_genome \
        --ref_index $ref_index \
        --ref_scaffold_name $ref_scaffold_name \
        --publish_dir $trim_align_output_dir
fi

if [ $call_variants = "yes" ]; then
    chkprevious "Step: call_variants" $trim_align_output_dir
    mkmissingdir $call_variants_output_dir
    nextflow -log ./.nextflow/nextflow.log \
        run ./pipeline/nextflow/call_variants.nf \
        -with-report $call_variants_output_dir/workflow_report.html \
        -profile $nextflow_profile \
        -resume \
        --cram_dir $trim_align_output_dir \
        --ref_genome $ref_genome \
        --ref_index $ref_index \
        --ref_scaffold_name $ref_scaffold_name \
        --ref_ploidy_file $ref_ploidy_file \
        --concatenate_vcf $concatenate_unfiltered_vcfs \
        --publish_dir $call_variants_output_dir
fi

if [ $filter_variants = "yes" ]; then
    chkprevious "Step: filter_variants" $call_variants_output_dir
    mkmissingdir $filter_variants_output_dir
    nextflow -log ./.nextflow/nextflow.log \
        run ./pipeline/nextflow/filter_variants.nf \
        -with-report $filter_variants_output_dir/workflow_report.html \
        -profile $nextflow_profile \
        -resume \
        --vcf_dir $call_variants_output_dir/chroms \
        --ref_index $ref_index \
        --ref_scaffold_name $ref_scaffold_name \
        --filtering_label $filtering_label \
        --min_alleles $filtering_min_alleles \
        --max_alleles $filtering_max_alleles \
        --max_missing $filtering_max_missing \
        --min_meanDP $filtering_min_meanDP \
        --max_meanDP $filtering_max_meanDP \
        --minDP $filtering_minDP \
        --maxDP $filtering_maxDP \
        --minQ $filtering_minQ \
        --keep $filtering_keep \
        --publish_dir $filter_variants_output_dir
fi

if [ $phase_variants = "yes" ]; then
    chkprevious "Step: phase_variants" $filter_variants_output_dir
    mkmissingdir $phase_variants_output_dir
    nextflow -log ./.nextflow/nextflow.log \
        run ./pipeline/nextflow/phase_variants.nf \
        -with-report $phase_variants_output_dir/workflow_report.html \
        -profile $nextflow_profile \
        -resume \
        --unphased_vcf ${filter_variants_output_dir}/variants_${filtering_label}.vcf.gz \
        --unphased_csi ${filter_variants_output_dir}/variants_${filtering_label}.vcf.gz.csi \
        --ref_index $ref_index \
        --ref_scaffold_name $ref_scaffold_name \
        --window_size $phasing_window_size \
        --ref_recombination_map_dir $ref_recombination_map_dir \
        --publish_dir $phase_variants_output_dir
fi

nextflow -log ./.nextflow/nextflow.log \
    run ./pipeline/nextflow/run_multiqc.nf \
    -profile $nextflow_profile \
    --results_dir $output_dir \
    --publish_dir $output_dir

# End work
