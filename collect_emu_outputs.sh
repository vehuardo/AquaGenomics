# Create the output file
output_file="concatenated_emu-outputs.tsv"

# Initialize a flag to track the first file
first_file=true

# Find and process each .tsv file
find . -type f -name "*_rel-abundance.tsv" | while read file; do
    # Extract the sample name
    sample_name=$(basename "$file" ._rel-abundance.tsv)

    if $first_file; then
        # Add header to the output file for the first file
        awk -v sample="$sample_name" 'NR==1 {print $0"\tSample"} NR>1 {print $0"\t"sample}' "$file" > "$output_file"
        first_file=false
    else
        # Append to the output file, skipping the header (NR>1)
        awk -v sample="$sample_name" 'NR>1 {print $0"\t"sample}' "$file" >> "$output_file"
    fi
done
