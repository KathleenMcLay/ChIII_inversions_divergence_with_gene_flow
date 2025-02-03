import os
import glob

# Define directory and output file
directory = "/QRISdata/Q6656/chapter_III/population_divergence"
output_file = "/QRISdata/Q6656/chapter_III/population_divergence/fst_all.txt"

# Find all .txt files in the directory
file_list = sorted(glob.glob(os.path.join(directory, "*weir.fst")))

if not file_list:
    print("No .txt files found in the directory.")
else:
    with open(output_file, "w") as outfile:
        first_file = True
        for file in file_list:
            with open(file, "r") as infile:
                lines = infile.readlines()
                if first_file:
                    outfile.writelines(lines)  # Write everything including header
                    first_file = False
                else:
                    outfile.writelines(lines[1:])  # Skip header after the first file
    print(f"Concatenation complete. Output saved to {output_file}")


import pandas as pd

# Define file paths
csv_file = "/QRISdata/Q6656/chapter_III/population_divergence/LD_all_data.csv"  # Update with your actual CSV file path
filter_file = "/home/uqkmcla4/ChIII_inversions_divergence_with_gene_flow/population_pair_divergence/popsLDfil.txt"  # Update with your actual text file path
output_file = "/QRISdata/Q6656/chapter_III/population_divergence/LD_all.csv"

# Load the filter list from the text file
with open(filter_file, "r") as f:
    filter_values = set(line.strip() for line in f if line.strip())

# Load the CSV file
df = pd.read_csv(csv_file)

# Filter the DataFrame
filtered_df = df[df['population'].isin(filter_values)]

# Save the filtered DataFrame to a new CSV file
filtered_df.to_csv(output_file, index=False)

print(f"Filtering complete. Output saved to {output_file}")