import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Function to create a list of barcode variations with one mismatch
def create_mismatch_variants(barcode):
    variants = [barcode]  # Include the perfect match
    for i in range(len(barcode)):
        for n in 'ACGT':
            if barcode[i] != n:
                variant = barcode[:i] + n + barcode[i+1:]
                variants.append(variant)
    return variants

# Function to trim barcode and additional 20 nucleotides from both ends of the sequence
def trim_barcode(record, barcode_seq, start_index, is_reverse=False, additional_trim=20):
    seq = str(record.seq)
    qual = record.letter_annotations["phred_quality"]
    
    if is_reverse:
        # If reverse barcode, trim sequence before the barcode and remove 20 nt from both ends
        trimmed_seq = Seq(seq[:start_index])  # Keep everything before the barcode
        trimmed_qual = qual[:start_index]
    else:
        # If forward barcode, trim sequence after the barcode and remove 20 nt from both ends
        trimmed_seq = Seq(seq[start_index + len(barcode_seq):])  # Keep everything after the barcode
        trimmed_qual = qual[start_index + len(barcode_seq):]

    # Now, remove 20 nucleotides from both ends of the remaining sequence
    if len(trimmed_seq) > 2 * additional_trim:  # Ensure there are enough bases left to trim
        final_seq = trimmed_seq[additional_trim:-additional_trim]
        final_qual = trimmed_qual[additional_trim:-additional_trim]
    else:
        final_seq = Seq("")  # Too short to trim, return an empty sequence
        final_qual = []

    return SeqRecord(
        final_seq,
        id=record.id,
        name=record.name,
        description=record.description,
        dbxrefs=record.dbxrefs,
        features=record.features,
        annotations=record.annotations,
        letter_annotations={"phred_quality": final_qual}
    )

# Function to find a barcode match allowing for one mismatch
def find_barcode(seq, barcode_variants):
    for variant in barcode_variants:
        start_index = seq.find(variant)
        if start_index != -1:
            return start_index  # Return the start index of the barcode
    return -1  # No barcode found

# Set up argument parser
parser = argparse.ArgumentParser(description='Demultiplex a FASTQ file based on custom barcodes.')
parser.add_argument('input_fastq', type=str, help='The FASTQ file to be demultiplexed')
args = parser.parse_args()

# Define barcodes and create mismatch variants
barcodes_fwd = {
    "BC01": "GGAGAAGAAGAA",
    "BC02": "AAGAGGCAAGAA",
    "BC03": "TCCTCCTAAGAA",
    "BC04": "CTCCGAAGAGAA",
    "BC05": "CCTATTGCAGAA",
    "BC06": "TGTCTGAAGGAA"
}
barcodes_rev = {
    "BC01": "TTCTTCTTCTCC",
    "BC02": "TTCTTGCCTCTT",
    "BC03": "TTCTTAGGAGGA",
    "BC04": "TTCTCTTCGGAG",
    "BC05": "TTCTGCAATAGG",
    "BC06": "TTCCTTCAGACA"    
}

# Generate mismatch variants
barcode_fwd_variants = {bc: create_mismatch_variants(seq) for bc, seq in barcodes_fwd.items()}
barcode_rev_variants = {bc: create_mismatch_variants(seq) for bc, seq in barcodes_rev.items()}

# Open output files for each barcode and unmatched reads
barcode_files = {bc: open(f"{bc}_output.fastq", "w") for bc in barcodes_fwd}
barcode_files["no_barcode"] = open("no-barcode.fastq", "w")

# Processing the FASTQ file
with open(args.input_fastq, "r") as fastq_file:
    for record in SeqIO.parse(fastq_file, "fastq"):
        seq = str(record.seq)
        matched = False

        # Check for forward barcodes within the first 100 nucleotides
        for bc_name, variants in barcode_fwd_variants.items():
            start_index = find_barcode(seq[:100], variants)
            if start_index != -1:
                print(f"Found forward barcode for {bc_name}")
                matched_record = trim_barcode(record, barcodes_fwd[bc_name], start_index)
                SeqIO.write(matched_record, barcode_files[bc_name], "fastq")
                matched = True
                break  # Stop looking if we've found a match

        # Only check the reverse if no forward barcode was found
        if not matched:
            for bc_name, variants in barcode_rev_variants.items():
                start_index = find_barcode(seq[-100:], variants)
                if start_index != -1:
                    print(f"Found reverse barcode for {bc_name}")
                    # For reverse barcodes, add the length of the sequence minus 100 to the start index
                    matched_record = trim_barcode(record, barcodes_rev[bc_name], len(seq) - 100 + start_index, is_reverse=True)
                    SeqIO.write(matched_record, barcode_files[bc_name], "fastq")
                    matched = True
                    break
        
        # Write unmatched reads to a separate file
        if not matched:
            SeqIO.write(record, barcode_files["no_barcode"], "fastq")

# Close output files
for bc, file in barcode_files.items():
    file.close()
