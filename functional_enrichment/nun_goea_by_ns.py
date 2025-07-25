# run_goea_by_ns.py
import sys
import os
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.associations import read_associations

def run_goea(study_file, pop_file, assoc_file, obo_file, label, output_dir):
    # Load data
    obodag = GODag(obo_file)

    # Read associations: geneid -> set of GO terms
    associations = read_associations(assoc_file, no_top=True)

    # Read population and study gene lists
    with open(pop_file) as f:
        pop_ids = [line.strip() for line in f if line.strip()]
    with open(study_file) as f:
        study_ids = [line.strip() for line in f if line.strip()]

    # GOEnrichmentStudyNS automatically splits by BP/MF/CC
    goea = GOEnrichmentStudyNS(
        pop_ids,
        associations,
        obodag,
        propagate_counts=False,
        methods=["fdr_bh"],
        alpha=0.05,
    )

    results = goea.run_study(study_ids)

    # Write results by namespace
    for ns in ['BP', 'MF', 'CC']:
        outfile = os.path.join(output_dir, f"go_enrichment_{label}_{ns}.tsv")
        goea.wr_tsv(outfile, results[ns])
        print(f"Wrote {ns} results to: {outfile}")

if __name__ == "__main__":
    study_file = sys.argv[1]
    pop_file = sys.argv[2]
    assoc_file = sys.argv[3]
    obo_file = sys.argv[4]
    label = sys.argv[5]
    output_dir = sys.argv[6]

    run_goea(study_file, pop_file, assoc_file, obo_file, label, output_dir)
