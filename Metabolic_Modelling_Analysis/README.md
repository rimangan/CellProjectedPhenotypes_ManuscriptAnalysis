# Metabolic modelling: Python scripts

**data_preprocessing.py** - Loads pseudobulk count matrices by cell type, CPM-normalizes, maps gene symbols to Ensembl.

**context_model.py** - Runs Riptide to contextualize Human-GEM with per-sample transcriptomes; prunes and saves one JSON model per sample.

**sampling.py** - Loads contextualized JSON models and runs OptGPSampler to generate flux samples per reaction; can split cell types into pickle files for parallel runs.

**sampling_analysis.py** - Loads median flux and metadata, runs rank-based regression with FDR and group comparisons (shuffle, CPP Tx groups, discordant, Path AD), and writes result CSVs per cell type.
