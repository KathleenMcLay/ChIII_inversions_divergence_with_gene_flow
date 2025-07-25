#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80GB
#SBATCH --job-name=Fun_enrichment
#SBATCH --time=3:00:00
#SBATCH --partition=general
#SBATCH --account=a_ortiz_barrientos_coe
#SBATCH --output=/home/uqkmcla4/scripts/Fun_enrichment.txt

# Load the SAMtools module (used for processing BAM files)
module load samtools

### Functional Enrichment Analysis ###

# Set base directory for storing intermediate and output files
BAS_DIR="/QRISdata/Q6656/chapter_III/functional_enrichment"

# Define input/output file paths
BAM="$BAS_DIR/mapped.bam"  # BAM file with transcript alignments
ANNOTATION="$BAS_DIR/annotated_without_contam_gene_ontology_terms.tsv" # GO annotation file

# Output files
POPULATION_FILE="$BAS_DIR/population.txt"       # All genes with GO annotations
ASSOCIATION_FILE="$BAS_DIR/association.txt"  # Gene-to-GO associations
OBO="/home/uqkmcla4/go-basic.obo"               # GO ontology file

### Step 1: Create the population and association files from annotation (all GO categories)

# awk '
# BEGIN {
#     FS = "\t"; OFS = "\t"
# }
# {
#     split($1, a, "|")
#     split(a[1], b, ".")
#     gene = b[1] "." b[2]  # Keep only first two parts
#     go = $2

#     if (gene != "" && go ~ /^GO:/) {
#         key = gene "_" go
#         if (!(key in seen)) {
#             go_terms[gene] = (gene in go_terms ? go_terms[gene] ";" go : go)
#             seen[key] = 1
#         }
#     }
# }
# END {
#     for (gene in go_terms) {
#         print gene > "'"$POPULATION_FILE"'"
#         print gene, go_terms[gene] > "'"$ASSOCIATION_FILE"'"
#     }
# }
# ' "$ANNOTATION"

### Step 2: Create per-inversion study files from named BEDs (all GO categories)

# for BED in "$BAS_DIR"/inversions_*.bed; do
#     BASENAME=$(basename "$BED")
#     NAME=$(echo "$BASENAME" | sed -E 's/^inversions_([^.]*)\.bed$/\1/')

#     GENE_ID_LIST="$BAS_DIR/gene_ids_${NAME}.txt"
#     MATCHED_GENES="$BAS_DIR/study_${NAME}.txt"

#     echo "Processing $BED -> $MATCHED_GENES"

#     samtools view -L "$BED" "$BAM" | \
#         awk '{print $1 "\t" $3 "\t" $4}' | \
#         sort -k2,2 -k3,3n | uniq > "$GENE_ID_LIST"

#     awk 'NR==FNR { ids[$1]; next }
#     {
#         split($1, parts, "|");
#         if (parts[3] in ids) {
#             split(parts[1], gene, ".");
#             print gene[1] "." gene[2];
#         }
#     }' "$GENE_ID_LIST" "$ANNOTATION" | sort -u > "$MATCHED_GENES"

#     echo "Done. Output written to $MATCHED_GENES"
# done

# ### Step 3: Create study files per inversion from coordinates in a single BED

# BED="$BAS_DIR/inversions.bed"  # BED file with inversion coordinates

# echo "Processing each inversion in $BED..."

# while IFS=$'\t' read -r scaffold start end; do
#     LABEL="${scaffold}_${start}_${end}"
#     TMP_BED="$BAS_DIR/tmp_${LABEL}.bed"
#     GENE_ID_LIST="$BAS_DIR/gene_ids_${LABEL}.txt"
#     MATCHED_GENES="$BAS_DIR/study_${LABEL}.txt"

#     echo -e "${scaffold}\t${start}\t${end}" > "$TMP_BED"

#     echo "Processing inversion $LABEL"

#     samtools view -L "$TMP_BED" "$BAM" | \
#         awk '{print $1 "\t" $3 "\t" $4}' | \
#         sort -k2,2 -k3,3n | uniq > "$GENE_ID_LIST"

#     awk 'NR==FNR { ids[$1]; next }
#     {
#         split($1, parts, "|");
#         if (parts[3] in ids) {
#             split(parts[1], gene, ".");
#             print gene[1] "." gene[2];
#         }
#     }' "$GENE_ID_LIST" "$ANNOTATION" | sort -u > "$MATCHED_GENES"

#     echo "Done: $MATCHED_GENES"

#     rm "$TMP_BED"
# done < "$BED"


### Step 4: Run GOATOOLS to identify enriched GO terms in inversions ###

#Run the GOATOOLS enrichment test

BAS_DIR="/QRISdata/Q6656/chapter_III/functional_enrichment"
POPULATION_FILE="$BAS_DIR/population.txt"
ASSOCIATION_FILE="$BAS_DIR/association_cleaned.txt"
OBO="/home/uqkmcla4/go-basic.obo"
OUTPUT_DIR="$BAS_DIR"

for STUDY in "$BAS_DIR"/study_D*.txt; do
    STUDY_BASENAME=$(basename "$STUDY")
    LABEL=$(echo "$STUDY_BASENAME" | sed -E 's/^study[._]?([^.]*)\.txt$/\1/')

    echo "Running GOATOOLS enrichment for $STUDY_BASENAME"

    python3 /home/uqkmcla4/scripts/run_goea_by_ns.py \
        "$STUDY" "$POPULATION_FILE" "$ASSOCIATION_FILE" "$OBO" "$LABEL" "$OUTPUT_DIR"
done        
