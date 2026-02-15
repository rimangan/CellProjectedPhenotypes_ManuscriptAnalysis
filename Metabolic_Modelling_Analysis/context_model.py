# Date: 15/02/2026
# This script is used to contextualize the models using the riptide model
# It is used to contextualize the models using the riptide model



import os
import pickle
import sys
import traceback
from functools import partial
from multiprocessing import Pool

import cobra
import numpy as np
import pandas as pd
from cobra.io import save_json_model

riptide_src_path = os.path.join(os.getcwd(), "riptide", "src")
if riptide_src_path not in sys.path:
    sys.path.insert(0, riptide_src_path)

# Import riptide
try:
    import riptide

    if not hasattr(riptide, "read_transcription_file"):
        from riptide import read_transcription_file, contextualize

        class RiptideModule:
            read_transcription_file = read_transcription_file
            contextualize = contextualize
        riptide = RiptideModule()
except ImportError:
    import importlib.util
    spec = importlib.util.spec_from_file_location("riptide", os.path.join(riptide_src_path, "riptide.py"))
    riptide = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(riptide)


if not hasattr(riptide, "read_transcription_file"):
    raise ImportError("Failed to import riptide.read_transcription_file. Check riptide installation.")

HumanGEM = cobra.io.read_sbml_model("Human-GEM.xml")

try:
    if "riptide" in sys.modules:
        import importlib
        importlib.reload(sys.modules["riptide"])
    else:
        import riptide
except Exception:
    import importlib.util
    riptide_path = os.path.join(riptide_src_path, "riptide.py")
    if os.path.exists(riptide_path):
        spec = importlib.util.spec_from_file_location("riptide", riptide_path)
        riptide = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(riptide)
    else:
        raise ImportError(f"Could not find riptide.py at {riptide_path}")

if hasattr(riptide, "read_transcription_file"):
    print("Riptide imported successfully (read_transcription_file available).")
else:
    print("riptide.read_transcription_file not found.")
    print(f"  Available attributes: {[x for x in dir(riptide) if not x.startswith('_')][:10]}")


def run_riptide_single_sample(row, cell_type, output_base_dir):
    sample_id = row["sample_id"]
    file_path = row["riptide_file"]

    safe_sample_id = (
        sample_id.replace("/", "_")
                 .replace("\\", "_")
                 .replace(" ", "_")
    )

    cell_out_dir = os.path.join(output_base_dir, cell_type)
    os.makedirs(cell_out_dir, exist_ok=True)

    try:
        print(f"Processing sample: {sample_id}")

        # Load transcriptome
        tpms = riptide.read_transcription_file(file_path)

        # Run riptide
        context_result = riptide.contextualize(
            model=HumanGEM,
            transcriptome=tpms,
            prune=True,
            fraction=0.8,
            objective=True,
            silent=False
        )

        model = context_result.model

        # Fix bounds
        for rxn in model.reactions:
            if hasattr(rxn.lower_bound, "evalf"):
                rxn.lower_bound = float(rxn.lower_bound.evalf())
            if hasattr(rxn.upper_bound, "evalf"):
                rxn.upper_bound = float(rxn.upper_bound.evalf())

        # Skip tiny models
        if len(model.reactions) < 50:
            return sample_id, {
                "warning": f"Model too small ({len(model.reactions)} reactions)",
                "riptide_file": file_path
            }

        # Save model
        model_outfile = os.path.join(
            cell_out_dir, f"{safe_sample_id}_model.json"
        )
        save_json_model(model, model_outfile)

        return sample_id, {
            "model_file": model_outfile,
            "flux_samples": context_result.flux_samples,
            "riptide_file": file_path
        }

    except Exception as e:
        return sample_id, {
            "error": str(e),
            "traceback": traceback.format_exc(),
            "riptide_file": file_path
        }


if __name__ == "__main__":

    OUTPUT_BASE_DIR = "riptide_output_biomass_objective"
    os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)

    cell_type = "Mic_P2RY12"

    manifest = pd.read_csv("riptide_manifest.csv")
    microglia_df = manifest[manifest["cell_type"] == f"{cell_type}"]

    if microglia_df.empty:
        raise ValueError(f"No {cell_type} samples found")

    riptide_results = {f"{cell_type}": {}}

    N_PROCESSES = 8
    print(f"Running with {N_PROCESSES} processes")

    worker = partial(
        run_riptide_single_sample,
        cell_type=f"{cell_type}",
        output_base_dir=OUTPUT_BASE_DIR
    )

    with Pool(processes=N_PROCESSES) as pool:
        for sample_id, result in pool.imap_unordered(
            worker,
            [row for _, row in microglia_df.iterrows()]
        ):
            riptide_results[f"{cell_type}"][sample_id] = result



    success = sum(
        1 for v in riptide_results[f"{cell_type}"].values()
        if "model_file" in v
    )

    print(f"{success}/{len(microglia_df)} samples completed successfully.")
    print(f"Riptide analysis completed for {cell_type}.")

