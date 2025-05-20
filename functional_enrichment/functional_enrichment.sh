#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=180GB
#SBATCH --job-name=Fun_enrichment
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --account=a_ortiz_barrientos_coe
#SBATCH --output=/home/uqkmcla4/scripts/Fun_enrichment.txt

# Load the SAMtools module (used for processing BAM files)
module load samtools

### Functional Enrichment Analysis ###

# Set base directory for storing intermediate and output files
BAS_DIR="/QRISdata/Q6656/chapter_III/functional_enrichment/"

### Step 1: Create the study file of genes that are present in the inversions ###

# Define input/output file paths for Step 1
BAM="$BAS_DIR/mapped.bam"                   # BAM file with transcript alignments
BED="$BAS_DIR/inversions.bed"              # BED file listing inversion regions
ANNOTATION="$BAS_DIR/annotated_without_contam_gene_ontology_terms.tsv" # GO annotation file
GENE_ID_LIST="$BAS_DIR/gene_ids.txt"       # Intermediate list of transcript IDs
MATCHED_GENES="$BAS_DIR/study.txt"         # Final study gene list (with biological_process GO terms)

# Extract reads from BAM that overlap inversion regions
echo "Extracting transcript IDs from BAM overlapping inversions..."
samtools view -L "$BED" "$BAM" | \
    awk '{print $1}' | sort | uniq > "$GENE_ID_LIST"

# Match the transcript IDs to GO annotations and extract only biological_process-related genes
echo "Matching GO annotations and extracting gene names..."
awk 'NR==FNR { ids[$1]; next } {
    split($1, parts, "|");
    if (parts[3] in ids && $4 == "biological_process") print parts[1];
}' "$GENE_ID_LIST" "$ANNOTATION" | sort | uniq > "$MATCHED_GENES"

# Notify completion of study file creation
echo "Done. Output written to $MATCHED_GENES"

### Step 2: Create the population and association files for all genes in the genome ###

# Input GO annotation file for the full genome
INPUT="$BAS_DIR/annotated_without_contam_gene_ontology_terms.tsv"

# Output files for GOATOOLS
POPULATION_FILE="$BAS_DIR/population.txt"      # List of all genes with biological_process GO terms
ASSOCIATION_FILE="$BAS_DIR/association.txt"    # Gene-to-GO term mappings (biological_process only)

# Parse the annotation file to build the population and association files
awk '
BEGIN {
    FS="\t"; OFS="\t"
}
# Only keep rows with the biological_process category
$4 == "biological Process" {
    split($1, a, "|")
    gene = a[1]
    go = $2
    if (gene != "" && go ~ /^GO:/)
        map[gene] = (map[gene] ? map[gene] ";" go : go)
}
# Write out the gene and associated GO terms
END {
    for (gene in map) {
        print gene > "'"$POPULATION_FILE"'"
        print gene, map[gene] > "'"$ASSOCIATION_FILE"'"
    }
}
' "$INPUT"

### Step 3: Run GOATOOLS to identify enriched GO terms in inversions ###

# Run the GOATOOLS enrichment test
# python3 /home/uqkmcla4/goatools/scripts/find_enrichment.py \
# "$BAS_DIR"/study.txt \                        # Study set: genes in inversions
# "$BAS_DIR"/population.txt \                  # Background gene set
# "$BAS_DIR"/association.txt \                 # GO term associations
# --obo /home/uqkmcla4/go-basic.obo \          # GO term definitions file
# --outfile "$BAS_DIR"/go_enrichment_results.xlsx \  # Output Excel file
# --pval 0.05 \                                 # P-value significance threshold
# --method fdr_bh                               # Multiple testing correction method (FDR Benjamini-Hochberg)
