module load samtools

### Functional Enrichment Analysis ###

BAS_DIR= "/QRISdata/Q6656/chapter_III/fun_enrichment/"

# 1_subset the list of gene_IDs & GO terms to just the ones that are present in the inversion

# convert .bam to .sam
samtools view -h "$BAS_DIR"/mapped.bam > "$BAS_DIR"/mapped.sam
# seperate the isoforms present in inversions
samtools view -h -b -L regions.bed "$BAS_DIR"/mapped.sam > "$BAS_DIR"/mapped_invonly.sam

#subset the list of go terms to just the ones that are present in the inversion
#extract the gene names from the output.sam file
awk '{print $3}' "$BAS_DIR"/mapped_invonly.sam | sort | uniq > "$BAS_DIR"/gene_names.txt

#extract the row from the .tsv file that contains the gene names
grep -f "$BAS_DIR"/gene_names.txt "$BAS_DIR"/annotated_without_contam_gene_ontology_terms.tsv > "$BAS_DIR"/go_terms_inversion.tsv

# Create a file that contains only unique gene names (study/population file for goatools)
awk '{print $1}' "$BAS_DIR"/go_terms_inversion.tsv | sort | uniq > "$BAS_DIR"/study.txt

#create a mapping of unique gene names to their associated GO terms. Output the final file with unique gene names in column 1 and associated GO terms in column 2 (association file for goatools)
awk 'BEGIN {FS="\t"; OFS="\t"} {go[$1] = (go[$1] ? go[$1] ";" $2 : $2)} END {for (gene in go) print gene, go[gene]}' "$BAS_DIR"/go_terms_inversion.tsv > "$BAS_DIR"/gene_to_go_terms.tsv
cat "$BAS_DIR"/gene_to_go_terms.tsv >> "$BAS_DIR"/association.txt

# 2_run GOAtools to identify GO terms that are enriched within inversions
# run GOAtools
python3 /home/uqkmcla4/goatools/scripts/find_enrichment.py data/study.txt data/population.txt data/association.txt --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --obo /home/uqkmcla4/go-basic.obo --outfile= out.xlsx 

