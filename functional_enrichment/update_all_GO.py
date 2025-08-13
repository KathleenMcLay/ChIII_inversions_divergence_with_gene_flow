import pandas as pd
import os
import glob
from goatools.obo_parser import GODag

# Configuration
input_directory = "/QRISdata/Q6656/chapter_III/5_gene_functional_enrichment/inversion_gene_files"  # Update this path
output_directory = "/QRISdata/Q6656/chapter_III/5_gene_functional_enrichment"  # Update this path
obo_file = "/home/uqkmcla4/go-basic.obo"

# Create output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Load GO DAG including obsolete terms
print("Loading GO DAG...")
go_dag = GODag(obo_file, 
               load_obsolete=True, 
               optional_attrs=['consider', 'replaced_by', 'alt_id'])

def map_current(go_id):
    term = go_dag.get(go_id)
    if term is None:
        return f"{go_id} (not_found)"
    if term.is_obsolete:
        if hasattr(term, 'replaced_by') and term.replaced_by:
            # Handle both list and set types
            replacement = list(term.replaced_by)[0] if isinstance(term.replaced_by, set) else term.replaced_by[0]
            return replacement
        elif hasattr(term, 'alt_id') and term.alt_id:
            # Handle both list and set types
            alt_id = list(term.alt_id)[0] if isinstance(term.alt_id, set) else term.alt_id[0]
            return alt_id
        elif hasattr(term, 'consider') and term.consider:
            # Handle both list and set types
            consider = list(term.consider)[0] if isinstance(term.consider, set) else term.consider[0]
            return consider
        else:
            return f"{go_id} (obsolete_no_replacement)"
    else:
        return go_id

def update_go_terms(go_column_str):
    if pd.isna(go_column_str):
        return ""
    terms = [tid.strip() for tid in go_column_str.split(',') if tid.strip()]
    updated = []
    
    for tid in terms:
        mapped_term = map_current(tid)
        # Only keep terms that don't have the obsolete_no_replacement marker
        if not mapped_term.endswith("(obsolete_no_replacement)"):
            updated.append(mapped_term)
    
    return ','.join(sorted(set(updated)))

def process_file(input_file_path, output_file_path):
    """Process a single TXT file"""
    try:
        print(f"Processing: {os.path.basename(input_file_path)}")
        
        # Read the file
        df = pd.read_csv(input_file_path, sep="\t")
        
        # Check if required columns exist
        required_columns = ["EggNOG.GO.Biological", "EggNOG.GO.Cellular", "EggNOG.GO.Molecular"]
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            print(f"  Warning: Missing columns {missing_columns} in {os.path.basename(input_file_path)}")
            print(f"  Available columns: {list(df.columns)}")
            return False
        
        # Update GO terms
        df["current.GO.Biological"] = df["EggNOG.GO.Biological"].apply(update_go_terms)
        df["current.GO.Cellular"]   = df["EggNOG.GO.Cellular"].apply(update_go_terms)
        df["current.GO.Molecular"]  = df["EggNOG.GO.Molecular"].apply(update_go_terms)
        
        # Save the updated file
        df.to_csv(output_file_path, sep="\t", index=False)
        print(f"  Successfully processed and saved to: {os.path.basename(output_file_path)}")
        return True
        
    except Exception as e:
        print(f"  Error processing {os.path.basename(input_file_path)}: {str(e)}")
        return False

# Main processing loop
def main():
    # Find all TXT files in the input directory
    txt_pattern = os.path.join(input_directory, "*.txt")
    txt_files = glob.glob(txt_pattern)
    
    if not txt_files:
        print(f"No TXT files found in {input_directory}")
        return
    
    print(f"Found {len(txt_files)} TXT files to process")
    
    processed_count = 0
    failed_count = 0
    
    for input_file in txt_files:
        # Generate output filename
        base_name = os.path.basename(input_file)
        name_without_ext = os.path.splitext(base_name)[0]
        output_file = os.path.join(output_directory, f"{name_without_ext}_with_current.txt")
        
        # Process the file
        if process_file(input_file, output_file):
            processed_count += 1
        else:
            failed_count += 1
    
    print(f"\nProcessing complete!")
    print(f"Successfully processed: {processed_count} files")
    print(f"Failed: {failed_count} files")

if __name__ == "__main__":
    main()