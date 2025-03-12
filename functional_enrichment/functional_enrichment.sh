module load samtools

### Functional Enrichment Analysis ###

BAS_DIR= "/QRISdata/Q6656/chapter_III/fun_enrichment/"

# per population pair (+coastal?) to capture each pops unique set of inversions. run goatools for all annotations 

# 1_subset the list of gene_IDs & GO terms to just the ones that are present in the inversion

# convert .bam to .sam
samtools view -h "$BAS_DIR"/mapped.bam > "$BAS_DIR"/mapped.sam

# create a list of unique gene names and associated GO terms for the whole genome (population.txt file for goatools)
awk 'BEGIN {FS="\t"; OFS="\t"} {go[$1] = (go[$1] ? go[$1] ";" $2 : $2)} END {for (gene in go) print gene, go[gene]}' "$BAS_DIR"/annotated_without_contam_gene_ontology_terms.tsv > "$BAS_DIR"/gene_to_go_terms.tsv

for file in ${directory}/sample_lists/*; do
    echo "current file is ${file}"
    filename=$(basename ${file} _samples.txt)

    # seperate the isoforms present in inversions
    samtools view -h -b -L ${filename}_regions.bed "$BAS_DIR"/mapped.sam > "$BAS_DIR"/${filename}_mapped_invonly.sam

    #subset the list of go terms to just the ones that are present in the inversion
    #extract the gene names from the output.sam file
    awk '{print $3}' "$BAS_DIR"/${filename}_mapped_invonly.sam | sort | uniq > "$BAS_DIR"/${filename}_gene_names.txt

    #extract the row from the .tsv file that contains the gene names
    grep -f "$BAS_DIR"/${filename}_gene_names.txt "$BAS_DIR"/annotated_without_contam_gene_ontology_terms.tsv > "$BAS_DIR"/${filename}_go_terms_inversion.tsv

    # Create a file that contains only unique gene names (study file for goatools) - study file is a list of genes that are present in the inversions ie. a subset of the full gene list across the whole genome
    awk '{print $1}' "$BAS_DIR"/${filename}_go_terms_inversion.tsv | sort | uniq > "$BAS_DIR"/${filename}_study.txt

    #create a mapping of unique gene names to their associated GO terms. Output the final file with unique gene names in column 1 and associated GO terms in column 2 (association file for goatools)
    awk 'BEGIN {FS="\t"; OFS="\t"} {go[$1] = (go[$1] ? go[$1] ";" $2 : $2)} END {for (gene in go) print gene, go[gene]}' "$BAS_DIR"/${filename}_go_terms_inversion.tsv > "$BAS_DIR"/${filename}_gene_to_go_terms.tsv
    cat "$BAS_DIR"/gene_to_go_terms.tsv >> "$BAS_DIR"/association.txt

    # 2_run GOAtools to identify GO terms that are enriched within inversions
    # run GOAtools
    python3 /home/uqkmcla4/goatools/scripts/find_enrichment.py data/study.txt data/population.txt data/association.txt --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --obo /home/uqkmcla4/go-basic.obo --outfile= ${filename}_goatools_out.xlsx 
done
