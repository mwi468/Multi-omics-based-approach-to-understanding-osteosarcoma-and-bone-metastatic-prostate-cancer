
###############################################################################################################
### adapter trimming in Parralel
#!/bin/bash

# Define the base directory containing the folders
base_dir="/mnt/e/ShanghaiOsteosarcoma_WES_RNAseq/HRA003260"
out_dir="/mnt/e/ShanghaiOsteosarcoma_WES_RNAseq/processed/trimmed"
log_file="/mnt/e/ShanghaiOsteosarcoma_WES_RNAseq/processed/trimmed/completed_samples.txt"
error_log="/mnt/e/ShanghaiOsteosarcoma_WES_RNAseq/processed/trimmed/error_samples.txt"

mkdir -p "$out_dir"

# Create the log files if they don't exist
touch "$log_file"
touch "$error_log"

## delete folders if they contain samples that were not completed 

mapfile -t completed_samples < "$log_file"

# Get a list of all folders in the output directory
for folder in "$out_dir"/*; do
  # Check if the folder is indeed a directory
  if [[ -d "$folder" ]]; then
    # Get the folder name (basename)
    folder_name=$(basename "$folder")
    
    # Check if the folder is in the list of completed samples
    if [[ ! " ${completed_samples[@]} " =~ " $folder_name " ]]; then
      # If not, delete the folder
      echo "Deleting folder: $folder"
      rm -rf "$folder"
    fi
  fi
done

# Function to process each sample
process_sample() {
    folder=$1
    folder_loc="${base_dir}/${folder}"

    # Define input FASTQ files for paired-end data with possible extensions
    f1_fq="${folder_loc}/${folder}_f1.fq.gz"      # Forward reads with _f1.fq.gz
    r2_fq="${folder_loc}/${folder}_r2.fq.gz"      # Reverse reads with _r2.fq.gz
    f1_fastq="${folder_loc}/${folder}_f1.fastq.gz"  # Forward reads with _f1.fastq.gz
    r2_fastq="${folder_loc}/${folder}_r2.fastq.gz"  # Reverse reads with _r2.fastq.gz

    # Check if the sample is already processed
    if grep -q "^${folder}$" "$log_file"; then
        echo "Sample ${folder} has already been processed. Skipping..."
        return
    fi

    # Check for _f1.fq.gz and _r2.fq.gz files, if not present check for _fastq.gz files
    if [ -f "$f1_fq" ] && [ -f "$r2_fq" ]; then
        echo "Files found for ${folder} with .fq.gz extension"
        f1="$f1_fq"
        r2="$r2_fq"
    elif [ -f "$f1_fastq" ] && [ -f "$r2_fastq" ]; then
        echo "Files found for ${folder} with .fastq.gz extension"
        f1="$f1_fastq"
        r2="$r2_fastq"
    else
        echo "Files not found for ${folder}" >> "$error_log"
        echo "Checked location: $folder_loc" >> "$error_log"
        echo "Files not found for ${folder}. Skipping..."
        return
    fi

    # Define output directory for this sample
    sample_output="${out_dir}/${folder}"
    mkdir -p "$sample_output"

    trimmed_f1="${sample_output}/${folder}_trimmed_f1.fq.gz"
    trimmed_r2="${sample_output}/${folder}_trimmed_r2.fq.gz"

    # Run Cutadapt to trim adapters
    cutadapt \
        -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
        -o "$trimmed_f1" -p "$trimmed_r2" \
        "$f1" "$r2" \
        -q 20,20 \
        --minimum-length 20 2>> "$error_log"

    # If the command succeeds, log the completed sample with a timestamp
    if [ $? -eq 0 ]; then
        echo "$(date +'%Y-%m-%d %H:%M:%S') - ${folder}" >> "$log_file"
        echo "Sample ${folder} processing completed successfully."
    else
        echo "$(date +'%Y-%m-%d %H:%M:%S') - Error processing ${folder}" >> "$error_log"
        echo "Error occurred while processing ${folder}. Check error log."
    fi
}

# Export the necessary variables and function for parallel execution
export base_dir out_dir log_file error_log
export -f process_sample

# List folders and run the script in parallel for 24 jobs at a time
folders=$(ls "$base_dir")

echo "$folders" | parallel -j 24 process_sample

#######################################################################################################
## multiqc checking in parralel

#!/bin/bash

# Define paths
input_file="RNAseqfiles2.tsv"
output_folder="/mnt/e/ShanghaiOsteosarcoma_WES_RNAseq/processed/fastqc"
num_jobs=20  # Number of parallel jobs

# Create the output folder if it doesn't exist
mkdir -p "$output_folder"

# Check if FastQC is installed
if ! command -v fastqc &> /dev/null; then
    echo "FastQC not found. Please install FastQC and ensure it is in your PATH."
    exit 1
fi

# Check if GNU Parallel is installed
if ! command -v parallel &> /dev/null; then
    echo "GNU Parallel not found. Please install GNU Parallel and ensure it is in your PATH."
    exit 1
fi

# Function to process each folder
process_folder() {
    local folder=$1
    local output_folder=$2

    # Get the base name of the folder to create specific output subfolders
    local sample_id=$(basename "$folder")

    # Define output directory for this sample
    local sample_output="$output_folder/$sample_id"
    mkdir -p "$sample_output"

    # Find FASTQ files using find
    local r1_files=$(find "$folder" -name '*f1.fq.gz' -print)
    local r2_files=$(find "$folder" -name '*r2.fq.gz' -print)

    # Debug output
    echo "Processing folder: $folder"
    echo "Files matching f1.fq.gz: $r1_files"
    echo "Files matching r2.fq.gz: $r2_files"

    # Check if any FASTQ files are found
    if [ -n "$r1_files" ] && [ -n "$r2_files" ]; then
        # Run FastQC on the FASTQ files
        fastqc $r1_files $r2_files -o "$sample_output"
        echo "Processed $sample_id and saved FastQC results in $sample_output"
    else
        echo "FASTQ files not found for $folder"
    fi
}

export -f process_folder

# Process folders in parallel
cat "$input_file" | parallel -j $num_jobs process_folder {} "$output_folder"

################################################################################################
### alignmet using STAR

#!/bin/bash

# Increase file limit size
ulimit -n 5000

# Define the directories and files
base_dir="/mnt/e/ShanghaiOsteosarcoma_WES_RNAseq/processed/trimmed"
out_dir="/mnt/e/ShanghaiOsteosarcoma_WES_RNAseq/processed/Aligned"
log_file="/mnt/e/ShanghaiOsteosarcoma_WES_RNAseq/processed/Aligned/completed_samples.txt"
error_log="/mnt/e/ShanghaiOsteosarcoma_WES_RNAseq/processed/Aligned/error_samples.txt"
genome_dir="/home/waleed/GenomeAnnotations"
tmp_dir_base="/mnt/ext4_partition/star_tmp"

# Create output and log directories if they don't exist
mkdir -p "$out_dir"
touch "$log_file"
touch "$error_log"

# Function to process each sample
process_sample() {
    folder=$1
    folder_loc="${base_dir}/${folder}"

    # Extract the 4th column from log_file and store unique completed samples in an array
    mapfile -t completed_samples < <(awk '{print $4}' "$log_file" | sort | uniq)

    # Check if the sample has already been processed by searching in the completed_samples array
    if [[ " ${completed_samples[@]} " =~ " ${folder} " ]]; then
        echo "Sample ${folder} has already been processed. Skipping..."
        return
    fi

    # Define input FASTQ files for paired-end data with possible extensions
    f1_fq="${folder_loc}/${folder}_trimmed_f1.fq.gz"
    r2_fq="${folder_loc}/${folder}_trimmed_r2.fq.gz"
    f1_fastq="${folder_loc}/${folder}_trimmed_f1.fastq.gz"
    r2_fastq="${folder_loc}/${folder}_trimmed_r2.fastq.gz"

    # Check for _f1.fq.gz and _r2.fq.gz files, if not present check for .fastq.gz files
    if [ -f "$f1_fq" ] && [ -f "$r2_fq" ]; then
        echo "Files found for ${folder} with .fq.gz extension"
        f1="$f1_fq"
        r2="$r2_fq"
    elif [ -f "$f1_fastq" ] && [ -f "$r2_fastq" ]; then
        echo "Files found for ${folder} with .fastq.gz extension"
        f1="$f1_fastq"
        r2="$r2_fastq"
    else
        echo "Files not found for ${folder}" >> "$error_log"
        echo "Checked location: $folder_loc" >> "$error_log"
        echo "Files not found for ${folder}. Skipping..."
        return
    fi

    # Define output directory for this sample
    sample_output="${out_dir}/${folder}"
    mkdir -p "$sample_output"

    # Create unique temporary directories for STAR
    tmp_dir1="${tmp_dir_base}/${folder}"
    rm -rf "$tmp_dir1"

    # Run STAR alignment
    STAR --genomeDir "$genome_dir" \
         --readFilesIn "$f1" "$r2" \
         --runThreadN 24 \
         --outFileNamePrefix "${sample_output}/" \
         --outSAMtype BAM SortedByCoordinate \
         --readFilesCommand zcat \
         --outTmpDir "$tmp_dir1"

    # Clean up tmp directory after run
    rm -rf "$tmp_dir1"

    # If the command succeeds, log the completed sample with a timestamp
    if [ $? -eq 0 ]; then
        echo "$(date +'%Y-%m-%d %H:%M:%S') - ${folder}" >> "$log_file"
        echo "Sample ${folder} processing completed successfully."
    else
        echo "$(date +'%Y-%m-%d %H:%M:%S') - Error processing ${folder}" >> "$error_log"
        echo "Error occurred while processing ${folder}. Check error log."
    fi
}

# Export the necessary variables and function for parallel execution
export base_dir out_dir log_file error_log genome_dir tmp_dir_base
export -f process_sample

# List folders and run the script in parallel for 1 job at a time
folders=$(ls "$base_dir")

echo "$folders" | parallel -j 1 process_sample

#######################################################################################################
## Feature Count Gene Expression 
#!/bin/bash

## usage first run ./featurecounts.sh

## usage second run so that function checks folder to skip samples already completed
# ./featurecounts10.sh --check-base /mnt/e/ShanghaiOsteosarcoma_WES_RNAseq/processed/Counts5

# Define input directory, where Aligned reads are in their sample folders
input_dir="/mnt/e/ShanghaiOsteosarcoma_WES_RNAseq/processed/Aligned"

# Define output directory base where output will be 
base_output_dir="/mnt/e/ShanghaiOsteosarcoma_WES_RNAseq/processed/Counts5"
output_dir="${base_output_dir}_$(date +%Y%m%d_%H%M%S)"  # Create a unique folder with a timestamp
log_file="$output_dir/completed_samples.log"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"
touch "$log_file"  # Create the log file if it doesn't exist

# Temporary directory for symbolic links with unique names
tmp_dir="$output_dir/tmp_symlinks"
mkdir -p "$tmp_dir"  # Create temp symlinks directory if it doesn't exist

# Initialize variables for checking base directory
check_base_dir=false
base_check_dir=""

# Parse arguments
if [[ "$1" == "--check-base" && -n "$2" ]]; then
    check_base_dir=true
    base_check_dir="$2"
    echo "Checking specified directory ($base_check_dir) for previously completed samples..."
elif [[ "$1" == "--check-base" && -z "$2" ]]; then
    echo "Error: Please specify a directory to check after '--check-base'."
    exit 1
fi

# Create a list of already completed samples for checking
mapfile -t completed_samples < <(awk '{print $1}' "$log_file" | sort | uniq)

# Create symlinks with unique names based on sample folder names
for bam in "$input_dir"/*/*.sortedByCoord.out.bam; do
    sample_name=$(basename "$(dirname "$bam")")
    symlink="$tmp_dir/${sample_name}_Aligned.sortedByCoord.out.bam"

    # Check if the sample is already in the completed_samples array or in the specified base directory
    if [[ " ${completed_samples[@]} " =~ " ${sample_name} " ]] || 
       ([[ "$check_base_dir" == true ]] && [ -e "${base_check_dir}/${sample_name}_counts.txt" ]); then
        echo "Counts for $sample_name already exist. Skipping..."
        continue
    fi

    # Create the symlink if it does not exist
    if [ ! -e "$symlink" ]; then
        ln -s "$bam" "$symlink"
    fi
done

# Function to run featureCounts on a single BAM file
run_featurecounts() {
    local bam_file=$1
    local output_dir=$2
    sample_name=$(basename "$bam_file" | cut -d '_' -f1)
    output_file="$output_dir/${sample_name}_counts.txt"

    # Skip if counts for this sample already exist in the log
    if grep -q "$sample_name" "$log_file"; then
        echo "Counts for $sample_name already exist. Skipping..."
        rm "$bam_file"  # Remove symlink immediately if already processed
        return
    fi

    # Run featureCounts for the current BAM file
    featureCounts -p --countReadPairs -a /home/waleed/GenomeAnnotations/gencode.v47.annotation.gtf \
    -o "$output_file" -T 2 -t exon -g gene_id "$bam_file"

    # Log completed sample to the log file
    echo "$sample_name" >> "$log_file"

    # Remove symlink after processing
    rm "$bam_file"
}

export -f run_featurecounts

# Batch processing to limit to 20 files at a time
bam_files=("$tmp_dir"/*.bam)
batch_size=20  # Limit to 20 files per batch

for ((i = 0; i < ${#bam_files[@]}; i += batch_size)); do
    batch=("${bam_files[@]:i:batch_size}")
    
    # Log the files in the current batch
    echo "Starting batch $((i / batch_size + 1)) with files:"
    for file in "${batch[@]}"; do
        echo " - $(basename "$file")"
    done
    
    # Process each batch in parallel with 10 jobs
    parallel -j 10 run_featurecounts ::: "${batch[@]}" ::: "$output_dir"
    
    # Confirm batch completion
    echo "Batch $((i / batch_size + 1)) completed."
done

# Clean up the temporary symlink directory
rm -r "$tmp_dir"
 

#########################################################################################
#Format files for combining
## only samples above 70% Mapping in feature counts were used to calculate TPM and combined

#!/bin/bash

# Define input and output directories
INPUT_DIR="/mnt/e/ShanghaiOsteosarcoma_WES_RNAseq/processed/Counts5/RNAsamplesfiltered"
OUTPUT_DIR="/mnt/e/ShanghaiOsteosarcoma_WES_RNAseq/TPM/FilteredColumns"
mkdir -p "$OUTPUT_DIR"

# Define the new header
NEW_HEADER="Geneid\tLength\tCount"

# Process each input file
for file in "$INPUT_DIR"/*_counts.txt; do
    sample_name=$(basename "$file" "_counts.txt")
    output_file="$OUTPUT_DIR/${sample_name}_filtered.txt"

    # Extract columns 1, 6, and 7, adding the header only once
    {
        echo -e "$NEW_HEADER"
        awk 'NR > 1 && $1 != "Geneid" {print $1, $6, $7}' OFS="\t" "$file"
    } > "$output_file"

    echo "Processed $file -> $output_file"
done

### Calculate TPM 

#!/bin/bash

# Input and output directories
input_dir="/mnt/e/ShanghaiOsteosarcoma_WES_RNAseq/TPM/FilteredColumns/"
output_dir="/mnt/e/ShanghaiOsteosarcoma_WES_RNAseq/TPM/TPMcalculated"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Function to calculate TPM for a single file
process_file() {
    input_file="$1"
    output_dir="$2"

    # Extract the base name without extension
    base_name=$(basename "$input_file")
    sample_name="${base_name%.*}"  # Remove extension
    output_file="$output_dir/${sample_name}_TPM.txt"

    # Declare an associative array for RPK values
    declare -A rpk_values
    total_rpk=0

    # Step 1: First pass to calculate RPK values
    while IFS=$'\t' read -r column1 column2 column3; do
        # Skip the header
        if [[ "$column1" == "Geneid" ]]; then
            continue
        fi

        # Store the second and third columns
        length="$column2"
        count="$column3"

        # Validate length and count values
        if [[ -z "$length" || -z "$count" ]]; then
            echo "Warning: Missing length or count for $column1 in $input_file" >&2
            continue
        fi

        # Calculate RPK (Reads Per Kilobase)
        if [[ "$length" -gt 0 ]]; then
            rpk=$(echo "$count / ($length / 1000)" | bc -l)
        else
            rpk=0
        fi

        # Store RPK value for each gene
        rpk_values["$column1"]=$rpk

        # Accumulate total RPK
        total_rpk=$(echo "$total_rpk + $rpk" | bc -l)
    done < "$input_file"

    # Step 2: Second pass to calculate TPM
    {
        echo -e "Geneid\tLength\tCount\tRPK\tTPM"  # Output header

        while IFS=$'\t' read -r column1 column2 column3; do
            # Skip the header
            if [[ "$column1" == "Geneid" ]]; then
                continue
            fi

            # Retrieve stored RPK value
            rpk="${rpk_values[$column1]}"

            # Calculate TPM
            if [[ "$total_rpk" != "0" ]]; then
                tpm=$(echo "($rpk / $total_rpk) * 1000000" | bc -l)
            else
                tpm=0
            fi

            # Output the results
            echo -e "$column1\t$column2\t$column3\t$rpk\t$tpm"
        done < "$input_file"
    } > "$output_file"

    echo "Processed $input_file -> $output_file"
}

export -f process_file

# Use GNU Parallel to process files in parallel
find "$input_dir" -type f | parallel -j "$(nproc)" process_file {} "$output_dir"

##### combine TPM files 

#!/bin/bash

# Input and output directories
input_dir="/mnt/e/ShanghaiOsteosarcoma_WES_RNAseq/TPM/TPMcalculated"
output_dir="/mnt/e/ShanghaiOsteosarcoma_WES_RNAseq/TPM/TPMsorted"
output_file2="${output_dir}/combined.txt"

mkdir -p "$output_dir"

# Process each _TPM.txt file in the input directory
for file in "$input_dir"/*_TPM.txt; do
    # Extract sample name and set output file
    sample_name=$(basename "$file" | cut -d'_' -f1)
    output_file="$output_dir/${sample_name}_TPM.txt"

    # Sort the file by Geneid (column 1) and extract only the TPM column (column 5), adding the header
    {
        echo -e "${sample_name}"  # Add the sample name as the header
        awk 'NR > 1 { print $1, $5 }' "$file" | sort -k1,1 | awk '{ print $2 }'
    } > "$output_file"

    echo "Processed $file -> $output_file"
done

paste -d '\t' "$output_dir"/*_TPM.txt >> "$output_file2"

first_file=$(ls "$input_dir"/*_TPM.txt | head -n 1)  # Get the first file
geneid_file="${output_dir}/Geneid.txt"           # Temporary file for Geneid

# Extract the Geneid column and write it to the temporary file
awk 'BEGIN {OFS="\t"} NR > 1 { print $1 }' "$first_file" | sort -k1,1 > "$geneid_file"
sed -i '1iGeneid' "$geneid_file" 

paste -d '\t' $geneid_file $output_file2 > "${output_dir}/final_combined.txt"

