#!/bin/bash
# Usage 
# After pipeline run:
# cd /home/pipeline/cidc_rnaseq
# ./rnaseq-schema-mapping.sh
# Wait for cp to complete and proceed with Portal upload

# Define a function to extract the ID from a filename or directory path
extract_id() {
    path=$1
    # Attempt to match the ID pattern (C followed by 8 alphanumeric characters)
    if [[ $path =~ C[A-Za-z0-9]{8} ]]; then
        id="${BASH_REMATCH[0]}"
    else
        id=""
    fi
    echo "$id"
}

# Define the source and target directory templates
declare -A file_map=(
    ["bam/{id}.bam"]="analysis/star/{id}.01/{id}.01.sorted.bam"
    ["bam/{id}.bam.bai"]="analysis/star/{id}.01/{id}.01.sorted.bam.bai"
    ["bam/{id}_sorted_bam_stats.txt"]="analysis/star/{id}.01/{id}.01.sorted.bam.stat.txt"
    ["bam/{id}.Aligned.toTranscriptome.out.bam"]="analysis/star/{id}.01/{id}.01.transcriptome.bam"
    ["bam/{id}.Chimeric.out.junction"]="analysis/star/{id}.01/{id}.01.Chimeric.out.junction"
    ["rseqc/read_distribution/{id}_read_distrib.txt"]="analysis/rseqc/read_distrib/{id}.01/{id}.01.txt"
    ["rseqc/tin_score/{id}_housekeeping.summary.txt"]="analysis/rseqc/tin_score/{id}.01/{id}.01.summary.txt"
    ["rseqc/tin_score/{id}_housekeeping.tin.xls"]="analysis/rseqc/tin_score/{id}.01/{id}.01.tin_score.txt"
    ["salmon/{id}/{id}_quant.sf"]="analysis/salmon/{id}.01/{id}.01.quant.sf"
    ["log/{id}_salmon.log"]="analysis/salmon/{id}.01/{id}.01.transcriptome.bam.log"
    ["salmon/{id}/aux_info/ambig_info.tsv"]="analysis/salmon/{id}.01/aux_info/ambig_info.tsv"
    ["salmon/{id}/aux_info/expected_bias.gz"]="analysis/salmon/{id}.01/aux_info/expected_bias.gz"
    ["salmon/{id}/aux_info/fld.gz"]="analysis/salmon/{id}.01/aux_info/fld.gz"
    ["salmon/{id}/aux_info/meta_info.json"]="analysis/salmon/{id}.01/aux_info/meta_info.json"
    ["salmon/{id}/aux_info/observed_bias.gz"]="analysis/salmon/{id}.01/aux_info/observed_bias.gz"
    ["salmon/{id}/aux_info/observed_bias_3p.gz"]="analysis/salmon/{id}.01/aux_info/observed_bias_3p.gz"
    ["salmon/{id}/cmd_info.json"]="analysis/salmon/{id}.01/cmd_info.json"
    ["salmon/{id}/logs/salmon_quant.log"]="analysis/salmon/{id}.01/logs/salmon_quant.log"
    ["centrifuge/{id}_centrifuge_addSample_report.tsv"]="analysis/microbiome/{id}.01/{id}.01_addSample_report.txt"
    ["trust4/{id}/{id}_report.tsv"]="analysis/trust4/{id}.01/{id}.01_report.tsv" 
  	["fusion/{id}/{id}.fusion_predictions.abridged_addSample.tsv"]="analysis/fusion/{id}.01/{id}.01.fusion_predictions.abridged_addSample.tsv"
    ["msisensor2/{id}_msisensor2.txt"]="analysis/msisensor/single/{id}.01/{id}.01_msisensor.txt"
    ["arcashla/{id}/{id}.genotype.json"]="analysis/neoantigen/{id}.01/{id}.01.genotype.json"
)

# Loop through the files and copy them to the target location
for src_template in "${!file_map[@]}"; do
    # Extract the destination template
    dest_template=${file_map[$src_template]}

    # Create a regex pattern for the source file, ensuring exact match with {id}
    src_pattern_regex=$(echo "$src_template" | sed 's/{id}/C[A-Za-z0-9]{8}/g')

    # Find all files that match the source template
    for src_file in $(find . -type f | grep -E "$src_pattern_regex$"); do
        # Extract the ID from the file or directory path
        id=$(extract_id "$src_file")

        if [ -z "$id" ]; then
            echo "No valid ID found for $src_file" >> error.log
            continue
        fi

        # Replace {id} in the destination path with the actual ID
        dest_file=${dest_template//"{id}"/"$id"}

        # Create the destination directory if it doesn't exist
        dest_dir=$(dirname "$dest_file")
        mkdir -p "$dest_dir"

        # Copy the file to the destination
        cp "$src_file" "$dest_file"

        echo "Copied $src_file to $dest_file"
    done
done