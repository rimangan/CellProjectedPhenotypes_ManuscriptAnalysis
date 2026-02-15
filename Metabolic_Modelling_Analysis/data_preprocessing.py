# Date: 15/02/2026
# This script is used to preprocess the data for the riptide model
# It is used to preprocess the data for the riptide model


import os
import pickle
import pandas as pd
import numpy as np


def normalize_cpm(expr_df, count_scaling_factor=1e6):
    """
    expr_df: samples × genes (raw counts)
    returns: samples × genes (CPM-normalized)
    """
    total_counts = expr_df.sum(axis=1)
    cpm = expr_df.div(total_counts, axis=0) * count_scaling_factor
    return cpm


BASE_INPUT_DIR = "../Metabolic_modelling/sum_pseudobulk_sub_cell_type"
BASE_OUTPUT_DIR = "riptide_input_files"

os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

with open("gene_dict.pkl", "rb") as f:
    gene_dict = pickle.load(f)

print(f"Loaded gene dictionary with {len(gene_dict)} mappings")


manifest_rows = []


cell_types = [
    d for d in os.listdir(BASE_INPUT_DIR)
    if os.path.isdir(os.path.join(BASE_INPUT_DIR, d))
]

print(f"Found {len(cell_types)} cell types")


for cell_type in cell_types:
    print(f"\nProcessing cell type: {cell_type}")
    cell_dir = os.path.join(BASE_INPUT_DIR, cell_type)

    # Output directory per cell type
    output_dir = os.path.join(BASE_OUTPUT_DIR, cell_type)
    os.makedirs(output_dir, exist_ok=True)

    # Identify metadata and count matrix
    csv_files = [f for f in os.listdir(cell_dir) if f.endswith(".csv")]

    metadata_file = next((f for f in csv_files if "meta" in f.lower()), None)
    count_file = next((f for f in csv_files if "count" in f.lower()), None)

    if metadata_file is None or count_file is None:
        print("Missing metadata or count matrix; skipping.")
        continue

    print(f"  Metadata: {metadata_file}")
    print(f"  Counts: {count_file}")

    metadata = pd.read_csv(os.path.join(cell_dir, metadata_file))
    sample_col_meta = metadata.columns[0]
    metadata[sample_col_meta] = metadata[sample_col_meta].astype(str)


    data = pd.read_csv(os.path.join(cell_dir, count_file))
    sample_col_counts = data.columns[0]
    data[sample_col_counts] = data[sample_col_counts].astype(str)

    # Keep only samples present in metadata
    valid_samples = set(metadata[sample_col_meta])
    data = data[data[sample_col_counts].isin(valid_samples)]

    print(f"  Samples retained: {data.shape[0]}")

    expr_df = data.set_index(sample_col_counts)
    expr_df = expr_df.astype(float)

    expr_cpm = normalize_cpm(expr_df)

    # Optional: remove genes with zero expression everywhere
    expr_cpm = expr_cpm.loc[:, expr_cpm.sum(axis=0) > 0]

    data_t = expr_cpm.T
    data_t.index.name = "gene_symbol"

    for sample_id in expr_cpm.index:
        sample_expression = data_t[sample_id]

        expression_df = pd.DataFrame({
            "gene_symbol": sample_expression.index,
            "expression": sample_expression.values
        })

        expression_df["ensembl_id"] = expression_df["gene_symbol"].map(gene_dict)
        expression_df = expression_df.dropna(subset=["ensembl_id"])

        riptide_df = expression_df[["ensembl_id", "expression"]]
        riptide_df.columns = ["gene_ID", "expression_value"]

        safe_sample_id = (
            sample_id.replace("/", "_")
                     .replace("\\", "_")
                     .replace(" ", "_")
        )

        output_file = os.path.join(
            output_dir, f"{safe_sample_id}_riptide_format.tsv"
        )

        riptide_df.to_csv(
            output_file,
            sep="\t",
            header=False,
            index=False
        )

        meta_row = metadata.loc[
            metadata[sample_col_meta] == sample_id
        ].iloc[0].to_dict()

        manifest_rows.append({
            "cell_type": cell_type,
            "sample_id": sample_id,
            "riptide_file": output_file,
            **meta_row
        })

    print(f"  Created {len(expr_cpm.index)} riptide files.")

manifest_df = pd.DataFrame(manifest_rows)
manifest_df.to_csv("riptide_manifest.csv", index=False)
