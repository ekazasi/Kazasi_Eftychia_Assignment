#!/bin/bash

# Function to load configuration values from a YAML file.
load_config() {
    local config_file="$1"

    if ! command -v yq &> /dev/null; then
        # Check if the 'yq' command is available; exit with an error if not.
        echo "Error: yq is required to parse YAML files. Install it and try again."
        exit 1
    fi

    # Extract configuration values from the YAML file using 'yq'.
    input_fasta=$(yq '.input_fasta' "$config_file")
    output_dir=$(yq '.output_directory' "$config_file")
    alignments_dir=$(yq '.alignments_directory' "$config_file")
    trees_dir=$(yq '.trees_directory' "$config_file")
    genes_dir=$(yq '.genes_directory' "$config_file")
    genes_trees_dir=$(yq '.genes_trees_directory' "$config_file")
    iqtree=$(yq '.iqtree' "$config_file") # Load iqtree binary path.
    iqtree_params=$(yq '.iqtree_params_first' "$config_file")
    iqtree_gene_params=$(yq '.iqtree_params_second' "$config_file")
    merge_script=$(yq '.merge_script' "$config_file")
    compare_script=$(yq '.compare_script' "$config_file")
    alignment_limit=$(yq '.alignment_limit' "$config_file")
    genes_limit=$(yq '.genes_limit' "$config_file")
    top5_exons=$(yq '.top5_exons_file' "$config_file")
    top5_genes=$(yq '.top5_genes_file' "$config_file")
    last5_exons=$(yq '.last5_exons_file' "$config_file")
    last5_genes=$(yq '.last5_genes_file' "$config_file")
}

# Function to extract alignments from the input FASTA file.
extract_alignments() {
    local input_fasta="$1"
    local output_dir="$2"
    local alignment_limit="$3"
    
    local count=0
    # Initialize a counter for the number of alignments extracted.

    local alignment=""
    # Variable to store the current alignment.

    while IFS= read -r line; do
        # Read each line from the input FASTA file.

        line=$(echo "$line" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        # Trim leading and trailing whitespace from the line.

        if [[ "$line" =~ ^">" ]]; then
            # If the line starts with '>', process it as a header.

            partial=$(echo "$line" | grep -o '_[^_]*_' | head -1 | sed 's/^_//;s/_$//')
            # Extract the partial identifier from the header.

            if [ -n "$partial" ]; then
                # If the partial identifier is not empty, update the header.
                line=">$partial"
            fi
        fi

        if [ -z "$line" ]; then
            # If the line is empty, save the current alignment to a file.

            if [ -n "$alignment" ]; then
                count=$((count+1))
                # Increment the alignment counter.

                if [ "$count" -gt "$alignment_limit" ]; then
                    # Stop processing if the alignment limit is reached.
                    echo "Reached limit of $alignment_limit alignments. Stopping."
                    break
                fi

                echo -e "$alignment" > "$output_dir/alignment_$count.fa"
                # Save the alignment to a file.

                alignment=""
                # Reset the alignment variable.
            fi
        else
            # If the line is not empty, append it to the current alignment.

            if [ -z "$alignment" ]; then
                alignment="$line"
            else
                alignment="$alignment\n$line"
            fi
        fi
    done < "$input_fasta"

    if [ -n "$alignment" ] && [ "$count" -lt "$alignment_limit" ]; then
        # Save the last alignment if it exists and the limit is not reached.
        count=$((count+1))
        echo -e "$alignment" > "$output_dir/alignment_$count.fa"
    fi

    echo "Extracted $((count - 1)) alignments to '$output_dir' directory."
}

# Function to extract unique genes from the input file.
extract_genes() {
    local input_file="$1"
    local output_dir="$2"
    local genes_limit="$3"
    

    local temp_dir="$output_dir/temp"
    mkdir -p "$temp_dir"
    # Create a temporary directory for intermediate files.

    declare -A unique_genes
    # Declare an associative array to store unique genes.

    gene_count=0
    # Initialize a counter for the number of unique genes.

    unique_genes_file="$temp_dir/unique_genes.txt" > "$unique_genes_file"
    # Create a file to store unique gene names.

    grep '^>' "$input_file" | while read -r header; do
        # Extract headers from the input file.

        gene=$(echo "$header" | cut -d'_' -f1 | tr -d '>')
        # Extract the gene name from the header.

        if [[ -z ${unique_genes[$gene]} ]]; then
            # If the gene is not already in the array, add it.

            unique_genes[$gene]=1
            echo "$gene" >> "$unique_genes_file"
            ((gene_count++))
        fi

        if (( gene_count >= genes_limit )); then
            # Stop processing if the gene limit is reached.
            break
        fi
    done

    output_file="$temp_dir/all_genes_chunk.fa" > "$output_file"
    # Create a file to store all gene chunks.

    while read -r gene; do
        # Process each unique gene.

        awk -v gene="$gene" '
        BEGIN { found = 0 }
        /^>/ {
            if ($0 ~ gene) {
                found = 1
                print $0
            } else {
                found = 0
            }
        }
        !/^>/ && found { print $0 }
        ' "$input_file" >> "$output_file"
        # Extract headers and sequences for the current gene.
    done < "$unique_genes_file"

    grouped_output_file="$temp/grouped_all_genes_chunk.fa" > "$grouped_output_file"
    # Create a file to store grouped gene chunks.

    while read -r header; do
        # Process each header in the output file.

        if [[ $header == ">"* ]]; then
            # If the header starts with '>', extract gene, species, and exon.

            gene=$(echo "$header" | cut -d'_' -f1 | tr -d '>')
            species=$(echo "$header" | cut -d'_' -f2)
            exon=$(echo "$header" | cut -d'_' -f3)
            read -r sequence
            # Read the sequence following the header.

            if [[ -n $sequence ]]; then
                # If the sequence is not empty, save it to the grouped file.
                echo -e "$gene\t$species\t$exon\t$sequence" >> "$grouped_output_file"
            fi
        fi
    done < "$output_file"

    while read -r line; do
        # Process each line in the grouped output file.

        gene=$(echo "$line" | cut -f1)
        # Extract the gene name from the line.

        echo -e "$line" >> "$output_dir/${gene}_grouped.fa"
        # Save the line to a file for the specific gene.
    done < "$grouped_output_file"

    for grouped_file in "$temp_dir"/*_grouped.fa; do
        # Process each grouped file in the temporary directory.

        declare -A species_sequences
        # Declare an associative array to store species sequences.

        species_order=()
        # Initialize an array to store the order of species.

        while read -r line; do
            # Process each line in the grouped file.

            species=$(echo "$line" | cut -f2)
            # Extract the species name from the line.

            sequence=$(echo "$line" | cut -f4)
            # Extract the sequence from the line.

            if [[ -z ${species_sequences["$species"]} ]]; then
                # If the species is not already in the array, add it.

                species_order+=("$species")
            fi

            species_sequences["$species"]+="$sequence"
            # Append the sequence to the species in the array.
        done < "$grouped_file"

        output_file="$output_dir/${grouped_file%.fa}_concatenated.fa" > "$output_file"
        # Create a file to store concatenated sequences.

        for species in "${species_order[@]}"; do
            # Process each species in the order array.

            echo -e "$species\n${species_sequences[$species]}" >> "$output_file"
            # Save the species and its concatenated sequence to the file.
        done

        unset species_sequences
        unset species_order
        # Clear the arrays for the next grouped file.
    done

    rm -rf "$temp_dir"
    # Remove the temporary directory.

    echo "Splitting complete. Files are saved in the '$split_output_dir' directory."
}

# Function to run IQ-TREE on a given file.
run_iqtree() {
    local file="$1"
    local output_dir="$2"
    local iqtree_params="$3"

    local filename=$(basename "$file")
    local output_prefix="$output_dir/${filename%.fa}"
    # Generate the output prefix for IQ-TREE results.

    echo "Running IQ-TREE on $file..."
    "$iqtree" -s "$file" $iqtree_params -pre "$output_prefix"
    # Run IQ-TREE with the specified parameters.

    if [ -f "${output_prefix}.treefile" ]; then
        echo "Tree file created: ${output_prefix}.treefile"
    else
        echo "Warning: Tree file not found for $file"
    fi
}

# Function to merge alignments and run IQ-TREE on the merged file.
merge_alignments() {
    local input_dir="$1"
    local output_folder="$2"
    local merge_script="$3"
    local iqtree_params="$4"
    

    local output_file="$input_dir/multi_alignment.fa"
    local tree_output_prefix="$output_folder/multi_alignment"
    # Define paths for the merged alignment file and tree output.

    if [ ! -d "$input_dir" ]; then
        # Check if the input directory exists.

        echo "Error: Directory $input_dir does not exist."
        exit 1
    fi

    ls "$input_dir"/*.fa | sort -V > "$input_dir/first_500.txt"
    # List and sort the first 500 alignment files.

    if ! command -v "$merge_script" &> /dev/null; then
        # Check if the merge script is available.

        echo "Error: merge_sequences script not found."
        exit 1
    fi

    "$merge_script" < "$input_dir/first_500.txt" > "$output_file"
    # Run the merge script to create the multi-sequence alignment.

    echo "Multi-sequence alignment saved to $output_file"

    echo "Running IQ-TREE on $output_file..."
    "$iqtree" -s "$output_file" $iqtree_params -pre "$tree_output_prefix"
    # Run IQ-TREE on the merged alignment file.

    if [ -f "${tree_output_prefix}.treefile" ]; then
        echo "Tree file created: ${tree_output_prefix}.treefile"
    else
        echo "Warning: Tree file not found for $output_file"
    fi

    rm -f "$input_dir/first_500.txt"
    # Remove the temporary file listing the first 500 alignments.
}

# Function to compare trees using Robinson-Foulds distances.
compare_trees() {
    local input_dir="$1"
    local compare_script="$2"
    local top5_file="$output_dir"/"$3"
    local last5_file="$output_dir"/"$4"
    

    local reference_tree="$input_dir/multi_alignment.treefile"
    local output_file="$output_dir/rf_distances.txt"
    # Define paths for the reference tree and output file.

    if [ ! -f "$reference_tree" ]; then
        # Check if the reference tree file exists.

        echo "Error: Reference tree '$reference_tree' not found."
        exit 1
    fi

    echo "Computing Robinson-Foulds distances..."
    > "$output_file"
    # Initialize the output file.

    for tree in "$input_dir"/*.treefile; do
        # Process each tree file in the input directory.

        if [[ "$tree" != "$reference_tree" ]]; then
            # Skip the reference tree file.

            Rscript "$compare_script" "$reference_tree" "$tree" >> "$output_file"
            # Compute the Robinson-Foulds distance and append to the output file.
        fi
    done

    sort -k2 -n "$output_file" | head -n 5 > "$top5_file"
    # Extract the top 5 closest trees.

    sort -k2 -n "$output_file" | tail -n 5 > "$last5_file"
    # Extract the 5 most different trees.

    rm "$output_dir/rf_distances.txt"
    # Remove the intermediate output file.

    echo "Top 5 closest trees saved in $top5_file"
    echo "The 5 most different trees saved in $last5_file"
}

 # Main function to execute the script's workflow.
main() {
   

    if [ "$#" -ne 1 ]; then
        # Ensure the script is called with exactly one argument (the config file).

        echo "Usage: $0 <config_file>"
        exit 1
    fi

    local config_file="$1"
    # Assign the first argument to the config_file variable.

    load_config "$config_file"
    # Load configuration values from the specified file.

    mkdir -p "$output_dir"
    mkdir -p "$alignments_dir" "$trees_dir" "$genes_dir" "$genes_trees_dir"
    
    # Working with alignments

    extract_alignments "$input_fasta" "$alignments_dir" "$alignment_limit"
    # Extract alignments from the input FASTA file.

    local fa_files=$(find "$alignments_dir" -name "*.fa" -type f)
    # Find all alignment files in the alignments directory.

    if [ -z "$fa_files" ]; then
        # Check if any alignment files were found.

        echo "Error: No .fa files found in $alignments_dir after extraction."
        exit 1
    fi

    for file in $fa_files; do
        # Run IQ-TREE on each alignment file.

        run_iqtree "$file" "$trees_dir" "$iqtree_params"
    done

    echo "All individual IQ-TREE analyses completed."

    merge_alignments "$alignments_dir" "$trees_dir" "$merge_script" "$iqtree_params"
    # Merge alignments and run IQ-TREE on the merged file.

    echo "All processes completed. Final tree from merged alignment is in: $trees_dir/multi_alignment.treefile"

    compare_trees "$trees_dir" "$compare_script" "$top5_exons" "$last5_exons"
    
    # Working with genes

    extract_alignments "$input_fasta" "$genes_dir" "$genes_limit"
    # Extract genes from the input FASTA file.

    local fa_files=$(find "$genes_dir" -name "*.fa" -type f)
    # Find all gene files in the genes directory.

    if [ -z "$fa_files" ]; then
        # Check if any gene files were found.

        echo "Error: No .fa files found in $genes_dir after extraction."
        exit 1
    fi

    for file in $fa_files; do
        # Run IQ-TREE on each gene file.

        run_iqtree "$file" "$genes_trees_dir" "$iqtree_params"
    done

    echo "All individual IQ-TREE analyses completed."

    merge_alignments "$genes_dir" "$genes_trees_dir" "$merge_script" "$iqtree_gene_params"
    # Merge genes and run IQ-TREE on the merged file.

    echo "All processes completed. Final tree from merged alignment is in: $genes_trees_dir/multi_alignment.treefile"

    compare_trees "$genes_trees_dir" "$compare_script" "$top5_genes" "$last5_genes"
    # Compare trees and save the results.
}

main "$@"
# Call the main function with all script arguments.
