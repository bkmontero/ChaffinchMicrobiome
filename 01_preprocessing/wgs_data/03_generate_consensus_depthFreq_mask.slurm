#!/bin/bash
#SBATCH -c 2
#SBATCH -n 1
#SBATCH --mem-per-cpu=8G
#SBATCH -t 24:00:00
#SBATCH --array=0-1  # 50 samples (0-49)



# Working directory
cd /mnt/netapp1/Store_IMIB/FringillaProject/Backup/nextflow_v2/consensus

# Input file list (49 cram files max)
INPUT_LIST=/mnt/netapp1/Store_IMIB/FringillaProject/Backup/nextflow_v2/modern_crams_1to50.list

# Reference genome
REF=/mnt/netapp1/Store_IMIB/FringillaProject/reference/bFriCoe1_complete.fa

# Ploidy file
PLOIDY=/mnt/netapp1/Store_IMIB/FringillaProject/reference/bFriCoe1.ploidy

# Output directory
FASTA_DIR=/mnt/netapp1/Store_IMIB/FringillaProject/Backup/nextflow_v2/consensus/fasta
mkdir -p $FASTA_DIR

# THRESHOLDS


MIN_DEPTH=3                    # Minimum depth (Pere uses 3)
MIN_QUAL=30                    # Minimum variant quality
MIN_ALLELE_FREQ=0.75           # Minimum allele frequency (Pere's 75% rule)

# SETUP

# Read input files
readarray -t INPUT_FILES < ${INPUT_LIST}
INPUT_FILE=${INPUT_FILES[${SLURM_ARRAY_TASK_ID}]}
INPUT_FILE=$(echo "${INPUT_FILE}" | tr -d '\r\n' | xargs)

if [[ -z "${INPUT_FILE}" || ! -f "${INPUT_FILE}" ]]; then
    echo "ERROR: Input file not found: ${INPUT_FILE}"
    exit 1
fi

# Extract sample name
SAMPLE=$(basename ${INPUT_FILE} | sed 's/\.[bc]ram$//' | sed 's/_[^_]*$//')

# Check file type
if [[ "${INPUT_FILE}" == *.cram ]]; then
    FILE_TYPE="CRAM"
else
    FILE_TYPE="BAM"
fi

echo "========================================="
echo "Creating consensus for: $SAMPLE"
echo "Input file: $INPUT_FILE"
echo "Min depth: ${MIN_DEPTH}X"
echo "Min allele frequency: ${MIN_ALLELE_FREQ}"
echo "Min quality: ${MIN_QUAL}"
echo "========================================="

# Create temporary directory
TMP_DIR=${SAMPLE}_tmp_${SLURM_JOBID}
mkdir -p $TMP_DIR
cd $TMP_DIR

# STEP 1: Generate mask (depth + allele frequency in one pass)

echo ""
echo "Step 1: Generating mask (depth + allele frequency filters)..."

samtools mpileup -f $REF -Q 0 -q 0 -aa $INPUT_FILE | \
    awk -v min_depth=$MIN_DEPTH \
        -v min_freq=$MIN_ALLELE_FREQ \
        'BEGIN {OFS="\t"}
        {
            chr = $1
            pos = $2
            ref = toupper($3)
            depth = $4
            bases = $5
            
            # Check depth first
            if (depth < min_depth) {
                print chr, pos-1, pos
                next
            }
            
            # If depth OK, check allele frequency
            if (depth > 0) {
                # Count bases
                count["A"] = 0
                count["C"] = 0
                count["G"] = 0
                count["T"] = 0
                
                # Clean pileup string
                bases_clean = bases
                
                # Remove read start markers (^Q)
                while (match(bases_clean, /\^./)) {
                    bases_clean = substr(bases_clean, 1, RSTART-1) substr(bases_clean, RSTART+2)
                }
                
                # Remove read end markers
                gsub(/\$/, "", bases_clean)
                
                # Convert to uppercase
                bases_clean = toupper(bases_clean)
                
                # Count bases character by character
                i = 1
                while (i <= length(bases_clean)) {
                    char = substr(bases_clean, i, 1)
                    
                    # Skip indels (+/-[0-9]+[ACGTN]+)
                    if (char == "+" || char == "-") {
                        i++
                        num_str = ""
                        while (i <= length(bases_clean) && match(substr(bases_clean, i, 1), /[0-9]/)) {
                            num_str = num_str substr(bases_clean, i, 1)
                            i++
                        }
                        indel_len = num_str + 0
                        i += indel_len
                        continue
                    }
                    
                    # Match to reference
                    if (char == "." || char == ",") {
                        count[ref]++
                    }
                    # Base call
                    else if (char == "A" || char == "C" || char == "G" || char == "T") {
                        count[char]++
                    }
                    
                    i++
                }
                
                # Find most common base
                max_count = 0
                for (base in count) {
                    if (count[base] > max_count) {
                        max_count = count[base]
                    }
                }
                
                # Calculate frequency
                total = count["A"] + count["C"] + count["G"] + count["T"]
                if (total > 0) {
                    freq = max_count / total
                    if (freq < min_freq) {
                        print chr, pos-1, pos
                    }
                }
            }
        }' > ${SAMPLE}_mask.bed

# Count masked sites
N_MASKED=$(wc -l < ${SAMPLE}_mask.bed)
GENOME_SIZE=$(awk '{sum+=$2} END {print sum}' ${REF}.fai)
PCT_MASKED=$(awk -v n=$N_MASKED -v g=$GENOME_SIZE 'BEGIN {printf "%.2f", (n/g)*100}')

echo "Sites masked: ${N_MASKED} (${PCT_MASKED}% of genome)"

# STEP 2: Call variants

echo ""
echo "Step 2: Calling variants..."

if [[ "$FILE_TYPE" == "CRAM" ]]; then
    bcftools mpileup \
        --fasta-ref $REF \
        -a DP,AD \
        -Ou \
        $INPUT_FILE | \
    bcftools call \
        --ploidy-file $PLOIDY \
        -mv \
        -O z \
        -o ${SAMPLE}_calls.vcf.gz
else
    bcftools mpileup \
        -a DP,AD \
        -Ou \
        -f $REF \
        $INPUT_FILE | \
    bcftools call \
        --ploidy-file $PLOIDY \
        -mv \
        -O z \
        -o ${SAMPLE}_calls.vcf.gz
fi

# STEP 3: Index VCF

echo ""
echo "Step 3: Indexing VCF..."
bcftools index ${SAMPLE}_calls.vcf.gz

# STEP 4: Normalize variants

echo ""
echo "Step 4: Normalizing variants..."
bcftools norm \
    -f $REF \
    ${SAMPLE}_calls.vcf.gz \
    -O z \
    -o ${SAMPLE}_calls.norm.vcf.gz

bcftools index ${SAMPLE}_calls.norm.vcf.gz

# STEP 5: Filter variants

echo ""
echo "Step 5: Filtering variants..."
echo "Filters: QUAL>${MIN_QUAL}, DP>=${MIN_DEPTH}, IndelGap=5bp"

bcftools filter \
    --IndelGap 5 \
    -i "QUAL>${MIN_QUAL} & FMT/DP>=${MIN_DEPTH}" \
    ${SAMPLE}_calls.norm.vcf.gz \
    -O z \
    -o ${SAMPLE}_filtered.vcf.gz

bcftools index ${SAMPLE}_filtered.vcf.gz

# Count variants
N_VARIANTS=$(bcftools view -H ${SAMPLE}_filtered.vcf.gz | wc -l)
echo "Variants passing filters: ${N_VARIANTS}"

# STEP 6: Create consensus with masking

echo ""
echo "Step 6: Creating consensus sequence..."

# Get sample name from VCF
SAMPLE_NAME=$(bcftools query -l ${SAMPLE}_filtered.vcf.gz | head -1)

# Create consensus with mask
bcftools consensus \
    -f $REF \
    -s $SAMPLE_NAME \
    -m ${SAMPLE}_mask.bed \
    ${SAMPLE}_filtered.vcf.gz | \
    bgzip -c > ${SAMPLE}.consensus.fa.gz

echo "Consensus sequence created"

# STEP 7: Index consensus

echo ""
echo "Step 7: Indexing consensus FASTA..."
samtools faidx ${SAMPLE}.consensus.fa.gz

# STEP 8: Extract mitochondrial DNA

echo ""
echo "Step 8: Extracting mitochondrial DNA..."

# Decompress for extraction
gunzip -c ${SAMPLE}.consensus.fa.gz > ${SAMPLE}.consensus.fa
samtools faidx ${SAMPLE}.consensus.fa

# Extract mtDNA
samtools faidx ${SAMPLE}.consensus.fa OY740764 > ${SAMPLE}_mtDNA.fa

# Rename header
sed -i "s/>OY740764/>${SAMPLE}_mtDNA/" ${SAMPLE}_mtDNA.fa

echo "Mitochondrial genome extracted"

# Remove uncompressed version
rm ${SAMPLE}.consensus.fa ${SAMPLE}.consensus.fa.fai

# STEP 9: Generate quality report

echo ""
echo "Step 9: Generating quality report..."

cat > ${SAMPLE}_consensus_report.txt << REPORT
CONSENSUS SEQUENCE QUALITY REPORT
Sample: ${SAMPLE}
Input file: ${INPUT_FILE}
File type: ${FILE_TYPE}
Date: $(date)

FILTERING PARAMETERS:
  Minimum depth: ${MIN_DEPTH}X
  Minimum allele frequency: ${MIN_ALLELE_FREQ}
  Minimum quality: ${MIN_QUAL}
  IndelGap filter: 5bp
  
RESULTS:
  Genome size: ${GENOME_SIZE} bp
  Sites masked: ${N_MASKED} (${PCT_MASKED}%)
  Variants called: ${N_VARIANTS}
  
MASKING METHOD:
  ✓ Sites with depth < ${MIN_DEPTH}X → N
  ✓ Sites where most common base < ${MIN_ALLELE_FREQ} → N
  ✓ Combined in single pass (no bedtools needed)
  
OUTPUT FILES:
  Full consensus: ${SAMPLE}.consensus.fa.gz
  Mitochondrial: ${SAMPLE}_mtDNA.fa
  Report: ${SAMPLE}_consensus_report.txt
REPORT

cat ${SAMPLE}_consensus_report.txt

# STEP 10: Copy to output and cleanup

echo ""
echo "Step 10: Copying files and cleaning up..."

# Copy essential files
cp ${SAMPLE}.consensus.fa.gz ${FASTA_DIR}/
cp ${SAMPLE}.consensus.fa.gz.fai ${FASTA_DIR}/
cp ${SAMPLE}_mtDNA.fa ${FASTA_DIR}/
cp ${SAMPLE}_consensus_report.txt ${FASTA_DIR}/

echo "Files copied to: ${FASTA_DIR}"

# Cleanup
TEMP_SIZE=$(du -sh . | cut -f1)
cd ..
rm -rf $TMP_DIR

echo "Temporary files removed (freed ~${TEMP_SIZE})"

# DONE
