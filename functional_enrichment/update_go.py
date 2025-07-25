from goatools.obo_parser import GODag

# === Config ===
obo_file = "/home/uqkmcla4/go-basic.obo"
assoc_file = "/QRISdata/Q6656/chapter_III/functional_enrichment/association.txt"
output_file = "/QRISdata/Q6656/chapter_III/functional_enrichment/association_cleaned.txt"
allow_consider = False  # Set to False to skip terms without replaced_by

# === Load GO Ontology DAG ===
print("Loading GO DAG...")
go_dag = GODag(obo_file, 
               load_obsolete=True, 
               optional_attrs=['consider', 'replaced_by', 'alt_id'])

# === Process association file ===
print("Processing:", assoc_file)

with open(assoc_file) as infile, open(output_file, "w") as outfile:
    for line in infile:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) != 2:
            continue
        gene, go_str = parts
        go_terms = go_str.split(";")

        updated_terms = set()

        for go_id in go_terms:
            print(f"Processing GO term: {go_id}")
            go_id = go_id.strip() # Remove any surrounding whitespace
            print(f"🔍 Checking GO term: {go_id}")
            if go_id not in go_dag:
                print(f"⚠️  Unknown GO term: {go_id} (skipped)")
                continue

            term = go_dag[go_id]
            if not term.is_obsolete:
                updated_terms.add(go_id)
            elif term.replaced_by:
                updated_terms.add(term.replaced_by[:10])
                print(f"🔁 Replaced {go_id} → {term.replaced_by[:10]}")
            elif allow_consider and term.consider:
                updated_terms.add(term.consider[:10])  # Use first consider term
                print(f"💡 Considered {go_id} → {term.consider[:10]}")
            else:
                print(f"❌ Skipped obsolete GO term with no replacement: {go_id}")

        if updated_terms:
            outfile.write(f"{gene}\t{';'.join(sorted(updated_terms))}\n")

print(f"\n✅ Cleaned association file written to: {output_file}")
