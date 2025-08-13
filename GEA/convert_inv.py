# convert inversion genotypes to bi-allelic counts for baypass

# import os
# import pandas as pd
# from collections import defaultdict

# input_dir = "/..."
# output_dir = "/..."

# # Make sure output directory exists
# os.makedirs(output_dir, exist_ok=True)

# for filename in os.listdir(input_dir):
#     if not filename.endswith(".txt"):
#         continue

#     filepath = os.path.join(input_dir, filename)
#     base = os.path.splitext(filename)[0]  # removes .txt
#     pop1 = base[:3]
#     pop2 = base[3:6]

#     df = pd.read_csv(filepath, sep="\t")

#     # Group rows by region
#     print(df.columns.tolist())
#     grouped = df.groupby("region")

#     results = []

#     for region, group in grouped:
#         counts = defaultdict(int)  # keys: ref_pop1, alt_pop1, ref_pop2, alt_pop2

#         for _, row in group.iterrows():
#             population = row["population"]
#             genotype = row["genotype"]

#             if pd.isnull(genotype):
#                 continue

#             try:
#                 genotype = int(genotype)
#             except ValueError:
#                 continue

#             if population == pop1:
#                 if genotype == 0:
#                     counts["ref_pop1"] += 2
#                 elif genotype == 1:
#                     counts["ref_pop1"] += 1
#                     counts["alt_pop1"] += 1
#                 elif genotype == 2:
#                     counts["alt_pop1"] += 2

#             elif population == pop2:
#                 if genotype == 0:
#                     counts["ref_pop2"] += 2
#                 elif genotype == 1:
#                     counts["ref_pop2"] += 1
#                     counts["alt_pop2"] += 1
#                 elif genotype == 2:
#                     counts["alt_pop2"] += 2

#         results.append([
#             counts["ref_pop1"],
#             counts["alt_pop1"],
#             counts["ref_pop2"],
#             counts["alt_pop2"]
#         ])

#     output_path = os.path.join(output_dir, f"{filename}_inversions.baypass")
#     with open(output_path, "w") as out:
#         for row in results:
#             out.write(" ".join(map(str, row)) + "\n")
#     print(f"Processed {filename} and saved to {output_path}")

import os
import pandas as pd
from collections import defaultdict

input_file = "/QRISdata/Q6656/chapter_III/inversion_genotypes_all.txt"
output_file = "/QRISdata/Q6656/chapter_III/inversion_bi-allelic_frequencies.txt"

# Define populations and column order
pop_list = ["D00", "H00", "D01", "H01", "D03", "H02", "D04", "H05", "D05", "H06", "D12", "H14"]

# Read input
df = pd.read_csv(input_file, sep="\t")

# Group by region
grouped = df.groupby("region")

# Construct header: region, then ref/alt per pop
header = ["region"] + [f"{pop}_R" for pop in pop_list] + [f"{pop}_A" for pop in pop_list]

# Container for output
output_rows = []

for region, group in grouped:
    counts = defaultdict(int)

    for _, row in group.iterrows():
        population = str(row["population"]).strip()
        genotype_raw = row["genotype"]

        if pd.isnull(genotype_raw) or population not in pop_list:
            continue

        try:
            genotype = int(float(genotype_raw))
        except ValueError:
            continue

        if genotype == 0:
            counts[f"ref_{population}"] += 2
        elif genotype == 1:
            counts[f"ref_{population}"] += 1
            counts[f"alt_{population}"] += 1
        elif genotype == 2:
            counts[f"alt_{population}"] += 2

    # Build output row
    row = [region]  # first column is region
    for pop in pop_list:
        row.append(counts[f"ref_{pop}"])
    for pop in pop_list:
        row.append(counts[f"alt_{pop}"])
    output_rows.append(row)

# Write output with headers
with open(output_file, "w") as out:
    out.write("\t".join(header) + "\n")
    for row in output_rows:
        out.write("\t".join(map(str, row)) + "\n")

print(f"Saved output to: {output_file}")
