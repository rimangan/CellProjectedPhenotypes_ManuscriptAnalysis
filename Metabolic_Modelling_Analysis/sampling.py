# 15/02/2026
# This script is used to sample the flux of the models using the OptGPSampler
# It is used to sample the flux of the models using the OptGPSampler


import os
import sys
import pickle
import pandas as pd
from cobra.io import load_json_model
from cobra.sampling import OptGPSampler
import os
import pandas as pd


RIPTIDE_SAMPLING_BASE = "riptide_output"
SPLIT_PICKLE_DIR = "."  # directory for part1/part2 pickle files

def get_all_cell_types():
    return sorted(
        d for d in os.listdir(RIPTIDE_SAMPLING_BASE)
        if os.path.isdir(os.path.join(RIPTIDE_SAMPLING_BASE, d))
    )

def split_and_save_cell_types():
    """Divide all cell types into two lists and save to pickle files."""
    all_ct = get_all_cell_types()
    n = len(all_ct)
    mid = (n + 1) // 2
    part1 = all_ct[:mid]
    part2 = all_ct[mid:]
    path1 = os.path.join(SPLIT_PICKLE_DIR, "cell_types_sampling_part1.pkl")
    path2 = os.path.join(SPLIT_PICKLE_DIR, "cell_types_sampling_part2.pkl")
    with open(path1, "wb") as f:
        pickle.dump(part1, f)
    with open(path2, "wb") as f:
        pickle.dump(part2, f)
    print(f"Saved {len(part1)} cell types to {path1}: {part1}")
    print(f"Saved {len(part2)} cell types to {path2}: {part2}")
    print("Run: python sampling.py cell_types_sampling_part1.pkl")
    print("Run: python sampling.py cell_types_sampling_part2.pkl")

if len(sys.argv) > 1:
    arg = sys.argv[1]
    if arg == "--split-and-save":
        split_and_save_cell_types()
        sys.exit(0)
    # treat as path to pickle file
    with open(arg, "rb") as f:
        cell_types = pickle.load(f)
    print(f"Loaded {len(cell_types)} cell types from {arg}: {cell_types}")
else:
    cell_types = get_all_cell_types()
    print(f"Found {len(cell_types)} cell types in {RIPTIDE_SAMPLING_BASE}: {cell_types}")


N_SAMPLES = 100
N_PROCESSES = 32


all_failed = []
for cell_type in cell_types:
    MODEL_BASE_DIR = os.path.join(RIPTIDE_SAMPLING_BASE, cell_type)
    OUTPUT_DIR = os.path.join("optgp_flux_samples", cell_type)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    flux_samples_dict = {}
    failed_samples = []

    model_files = sorted(
        f for f in os.listdir(MODEL_BASE_DIR)
        if f.endswith("_model.json")
    )

    print("\n" + "=" * 80)
    print(f"Running OptGPSampler for {cell_type} models ({len(model_files)} models)")
    print("=" * 80)

    for idx, model_file in enumerate(model_files, 1):
        sample_id = model_file.replace("_model.json", "")
        model_path = os.path.join(MODEL_BASE_DIR, model_file)
       

        out_file = os.path.join(OUTPUT_DIR, f"{sample_id}_optgp_flux_samples.pkl")

        print(f"\n[{idx}/{len(model_files)}] {cell_type} | Sample: {sample_id}")

        if os.path.exists(out_file):
            print(f" Already sampled, skipping")
            continue

        try:
            model = load_json_model(model_path)
            sampler = OptGPSampler(model, processes=N_PROCESSES)
            samples = sampler.sample(N_SAMPLES)
            samples.to_pickle(out_file)
            flux_samples_dict[sample_id] = samples
            print(f"  Sampled {samples.shape[0]} points, {samples.shape[1]} reactions → {out_file}")
        except Exception as e:
            print(f"Failed sampling {sample_id}: {str(e)}")
            failed_samples.append(sample_id)
            all_failed.append((cell_type, sample_id))

    if flux_samples_dict:
        summary_path = os.path.join(OUTPUT_DIR, f"optgp_flux_samples_{cell_type}.pkl")
        with open(summary_path, "wb") as f:
            pickle.dump(flux_samples_dict, f)

    print(f"\n{cell_type}: newly sampled={len(flux_samples_dict)}, failed={len(failed_samples)}, skipped={len(model_files) - len(flux_samples_dict) - len(failed_samples)}")

print("\n" + "=" * 80)
print("OptGPSampler finished for all cell types")
print(f"Total failed across all cell types: {len(all_failed)}")


