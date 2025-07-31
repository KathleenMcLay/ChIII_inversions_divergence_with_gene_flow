# convert inversion genotypes to bi-allelic counts for baypass

import os
import pandas as pd
from collections import defaultdict

input_dir = "/..."
output_dir = "/..."

# Make sure output directory exists
os.makedirs(output_dir, exist_ok=True)

for filename in os.listdir(input_dir):
    if not filename.endswith(".txt"):
        continue

    filepath = os.path.join(input_dir, filename)
    base = os.path.splitext(filename)[0]  # removes .txt
    pop1 = base[:3]
    pop2 = base[3:6]

    df = pd.read_csv(filepath, sep="\t")

    # Group rows by region
    print(df.columns.tolist())
    grouped = df.groupby("region")

    results = []

    for region, group in grouped:
        counts = defaultdict(int)  # keys: ref_pop1, alt_pop1, ref_pop2, alt_pop2

        for _, row in group.iterrows():
            population = row["population"]
            genotype = row["genotype"]

            if pd.isnull(genotype):
                continue

            try:
                genotype = int(genotype)
            except ValueError:
                continue

            if population == pop1:
                if genotype == 0:
                    counts["ref_pop1"] += 2
                elif genotype == 1:
                    counts["ref_pop1"] += 1
                    counts["alt_pop1"] += 1
                elif genotype == 2:
                    counts["alt_pop1"] += 2

            elif population == pop2:
                if genotype == 0:
                    counts["ref_pop2"] += 2
                elif genotype == 1:
                    counts["ref_pop2"] += 1
                    counts["alt_pop2"] += 1
                elif genotype == 2:
                    counts["alt_pop2"] += 2

        results.append([
            counts["ref_pop1"],
            counts["alt_pop1"],
            counts["ref_pop2"],
            counts["alt_pop2"]
        ])

    output_path = os.path.join(output_dir, f"{filename}_inversions.baypass")
    with open(output_path, "w") as out:
        for row in results:
            out.write(" ".join(map(str, row)) + "\n")
    print(f"Processed {filename} and saved to {output_path}")

### for a single input file with all populations ###

import os
import pandas as pd
from collections import defaultdict

input_file = "/.../inversion_genotypes.txt"  # replace with actual file path
output_dir = "/.../GEA_inversions_all"
os.makedirs(output_dir, exist_ok=True)

# Define populations and output column order
pop_list = ["D00", "H00", "D01", "H01", "D03", "H02", "D04", "H05", "D12", "H14"]

# Read the input file
df = pd.read_csv(input_file, sep="\t")

# Group by region
grouped = df.groupby("region")

# Container for output
output_rows = []

for region, group in grouped:  # for each inversion
    print(region)
    counts = defaultdict(int)  # keys like ref_D00, alt_D00, ...

    for _, row in group.iterrows():
        population = row["population"]
        print(population)
        genotype = row["genotype"]
        print(genotype)
        if pd.isnull(genotype) or population not in pop_list:
            print("no population match")
            continue
        try:
            genotype = int(genotype)
        except ValueError:
            continue

        if genotype == 0:
            counts[f"ref_{population}"] += 2
        elif genotype == 1:
            counts[f"ref_{population}"] += 1
            counts[f"alt_{population}"] += 1
        elif genotype == 2:
            counts[f"alt_{population}"] += 2

    # Create output row with region name first
    row = [region]  # add region name
    for pop in pop_list:
        row.append(counts[f"ref_{pop}"])
        row.append(counts[f"alt_{pop}"])

    output_rows.append(row)

# Write final output file
output_path = os.path.join(output_dir, "all_pops_inversions.baypass")
with open(output_path, "w") as out:
    for row in output_rows:
        out.write(" ".join(map(str, row)) + "\n")

print(f"Saved BayPass-formatted output to {output_path}")
