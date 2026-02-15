import os
import pickle
import pandas as pd
import statsmodels.api as sm
from scipy.stats import rankdata, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from cobra.io import read_sbml_model

model = read_sbml_model("../Human-GEM.xml")



def get_reaction_info(rxn_id):
    try:
        reaction = model.reactions.get_by_id(rxn_id)
        
        # Compartments
        compartments = {met.compartment for met in reaction.metabolites}
        compartments_names = [model.compartments[c] for c in compartments]
        compartments_str = ", ".join(compartments_names)
        
        # Reaction type
        reaction_type = "Metabolic" if len(compartments) == 1 else "Transport"
        
        # Subsystem
        subsystem = reaction.subsystem if hasattr(reaction, 'subsystem') else "NA"
        
        # Metabolites involved (use names instead of IDs)
        metabolites = [met.name if hasattr(met, 'name') and met.name else met.id for met in reaction.metabolites]
        metabolites_str = ", ".join(metabolites) if metabolites else "NA"
        
        return pd.Series([compartments_str, reaction_type, subsystem, metabolites_str])
    
    except KeyError:
        return pd.Series(["NA", "Unknown", "NA", "NA"])


def rank_regression_test(
    merged_metadata,
    sampling_aggregated_df,
    
    group1,
    group2,
    group_col = None,
    covariates=None,
    categorical_covariates=None,
    id_col="projid"
):
    """
    Rank-based regression to test group effect on each reaction,
    adjusting for covariates (categorical or continuous).

    Parameters:
        merged_metadata : pd.DataFrame
            Metadata with group and covariate information
        sampling_aggregated_df : pd.DataFrame
            Reactions x samples (rows must match merged_metadata)
        group_col : str
            Column name for group
        group1, group2 : str
            Names of the two groups to compare
        covariates : list of str
            Continuous covariates to include (optional)
        categorical_covariates : list of str
            Categorical covariates to include (optional)
        id_col : str
            Column identifying samples (optional)

    Returns:
        pd.DataFrame : reaction-level statistics
    """

    results = []

    # --- 1. Filter metadata for the two groups ---
    meta = merged_metadata[merged_metadata[group_col].isin([group1, group2])].copy()

    # Binary group variable: 0 = group1, 1 = group2
    meta["group_binary"] = (meta[group_col] == group2).astype(int)

    # --- 2. Dummy encode categorical covariates ---
    if categorical_covariates is not None and len(categorical_covariates) > 0:
        meta = pd.get_dummies(meta, columns=categorical_covariates, drop_first=True)

        # Update covariates list to include new dummy columns
        if covariates is None:
            covariates = []
        dummy_cols = [c for c in meta.columns if any(c.startswith(cat + "_") for cat in categorical_covariates)]
        covariates += dummy_cols

    # --- 3. Loop over reactions ---
    for rxn in sampling_aggregated_df.columns:

        y = sampling_aggregated_df.loc[meta.index, rxn]

        if y.notna().sum() < 10:
            continue

        y_rank = rankdata(y)

        # Build design matrix
        X_cols = ["group_binary"] + (covariates if covariates is not None else [])
        X = meta[X_cols]
        X = X.astype(float)
     
        X = sm.add_constant(X)

        # Fit OLS on ranks
        model = sm.OLS(y_rank, X).fit(cov_type="HC3")  # robust SE

        results.append({
            "reaction": rxn,
            "coef_group": model.params["group_binary"],
            "se_group": model.bse["group_binary"],
            "t_group": model.tvalues["group_binary"],
            "p_value": model.pvalues["group_binary"],
            "n": len(y)
        })

    results_df = pd.DataFrame(results)

    return results_df




def compare_groups_whitney(
    merged_metadata,
    sampling_aggregated_df,
    group1,
    group2,
    id_col="projid",
    group_col="subgroup",
    alternative="two-sided"
):
    """
    Compare two subgroups using Mann–Whitney U test on flux medians.

    Returns:
        results_df : reaction-level statistics
        group1_df  : filtered median matrix for group1
        group2_df  : filtered median matrix for group2
    """

    # --- 1. Select individuals ---
    group_1_metadata = merged_metadata[merged_metadata[f'{group_col}'].isin([group1])]
    group_2_metadata = merged_metadata[merged_metadata[f'{group_col}'].isin([group2])]
    group1_df = sampling_aggregated_df.loc[group_1_metadata.index]
    group2_df = sampling_aggregated_df.loc[group_2_metadata.index]
    results = []
    for rxn in sampling_aggregated_df.columns:
        x = group1_df[rxn]
        y = group2_df[rxn]

        # Skip if not enough data
        if len(x) < 5 or len(y) < 5:
            continue

        stat, p = mannwhitneyu(x, y, alternative=alternative)

        results.append({
            "reaction": rxn,
            "group1": group1,
            "group2": group2,
            "n_group1": len(x),
            "n_group2": len(y),
            "median_group1": x.median(),
            "median_group2": y.median(),
            "effect_median_diff": y.median() - x.median(),
            "p_value": p
        })

    results_df = pd.DataFrame(results)

    return results_df, group1_df, group2_df


cell_done = []
cell_error = []
remove_cell_type = []

for cell_type in os.listdir("../new_opt_flux/opt_flux/"):
    if cell_type in cell_done:
        print(f"Cell type already done: {cell_type}")
        continue

    if cell_type in remove_cell_type:
        print(f"Skipping cell type: {cell_type}")
        continue

    meta_data = pd.read_csv("merged_df_final_12_18.csv", index_col=0)
    
    meta_data_2 = pd.read_csv("meta_data_with_batch.csv")
    
    merged_metadata = meta_data.merge(
        meta_data_2,
        on="projid",
        how="inner"
    )
    batch_counts = merged_metadata["batch_x"].value_counts()
    bottom_batches = batch_counts.nsmallest(5).index.tolist()
    merged_metadata = merged_metadata[~merged_metadata["batch_x"].isin(bottom_batches)].copy()
    samples = {}
    folder = f"../new_opt_flux/optgp_flux_samples/{cell_type}"
    for fname in os.listdir(folder):
        if fname == f"optgp_flux_samples_{cell_type}.pkl":
            continue
        if fname.endswith((".pkl", ".pickle")):
            key = os.path.splitext(fname)[0]
            with open(os.path.join(folder, fname), "rb") as f:
                    samples[key.split("_")[0]] = pickle.load(f)
     
    
    median_flux = {}
    
    for sample_id, sampling in samples.items():
        # median per reaction
        median_flux[sample_id] = sampling.median(axis=0)
    
    # Combine into one DataFrame
    median_flux_df = pd.DataFrame(median_flux).T
    
    median_flux_df = median_flux_df.fillna(0)

    mean_flux = {}
    
    
    keep_subject = list(set(median_flux_df.index).intersection(set(merged_metadata.subject)))
    median_flux_df = median_flux_df.loc[keep_subject]
    merged_metadata = merged_metadata.set_index('subject')
    merged_metadata  = merged_metadata.loc[median_flux_df.index]
    
    
    numerical_cols = merged_metadata.select_dtypes(include=["int64", "float64"]).columns.tolist()
    categorical_cols = merged_metadata.select_dtypes(include=["object", "category", "bool"]).columns.tolist()
    merged_metadata["msex"] = merged_metadata["msex"].astype("category")

    try:
        results_df_shuffle_group= rank_regression_test(
            merged_metadata,
            sampling_aggregated_df = median_flux_df,
            group1="GroupA",
            group2="GroupB",
            group_col = 'shuf_group',
            covariates= ['pmi', 'age_death'],
            categorical_covariates=['batch_x', 'msex'],
            id_col="projid")
        
        results_df_shuffle_group["q_value"] = multipletests(
            results_df_shuffle_group["p_value"],
            method="fdr_bh"
        )[1] 
        
        results_df_shuffle_group = results_df_shuffle_group.sort_values(by ='q_value')
        results_df_shuffle_group["reaction_name"] = results_df_shuffle_group["reaction"].map(
            lambda x: model.reactions.get_by_id(x).name
            )
        results_df_shuffle_group[["compartments_str", "reaction_type", "subsystem", "metabolites"]] = results_df_shuffle_group["reaction"].apply(get_reaction_info)

        base_dir = "results_new"
        output_dir = os.path.join(base_dir, cell_type, "shuffle_group")
        os.makedirs(output_dir, exist_ok=True)
        
        results_df_shuffle_group.to_csv(
            os.path.join(output_dir, "rank_regression_fdr05.csv"),
            index=False
        )
        results_df_shuffle_group = results_df_shuffle_group[results_df_shuffle_group.q_value <= 0.05] 
    
        print(f"Shuffle: {len(results_df_shuffle_group)} reactions at FDR 0.05")

        results_CPP_based_Tx_subgroups_regression = rank_regression_test(
            merged_metadata,
            sampling_aggregated_df = median_flux_df,
            group1="PathNonAD-TxNonAD",
            group2="PathAD-TxAD",
            group_col = 'subgroup',
            covariates= ['pmi', 'age_death'],
            categorical_covariates=['batch_x','msex'],
            id_col="projid")
        
        results_CPP_based_Tx_subgroups_regression["q_value"] = multipletests(
            results_CPP_based_Tx_subgroups_regression["p_value"],
            method="fdr_bh"
        )[1]
        
        results_CPP_based_Tx_subgroups_regression = results_CPP_based_Tx_subgroups_regression.sort_values(by ='q_value')
        results_CPP_based_Tx_subgroups_regression["reaction_name"] = results_CPP_based_Tx_subgroups_regression["reaction"].map(
            lambda x: model.reactions.get_by_id(x).name
        )
        results_CPP_based_Tx_subgroups_regression
        results_CPP_based_Tx_subgroups_regression[["compartments_str", "reaction_type", "subsystem", "metabolites"]] = results_CPP_based_Tx_subgroups_regression["reaction"].apply(get_reaction_info)

        output_dir = os.path.join(base_dir, cell_type, "cpp_based_tx_group")
        os.makedirs(output_dir, exist_ok=True)
        
        results_CPP_based_Tx_subgroups_regression.to_csv(
            os.path.join(output_dir, "rank_regression_fdr05.csv"),
            index=False
        )
        
        results_CPP_based_Tx_subgroups_regression = results_CPP_based_Tx_subgroups_regression[results_CPP_based_Tx_subgroups_regression.q_value <= 0.05]
        print(f"CPP: {len(results_CPP_based_Tx_subgroups_regression)} reactions at FDR 0.05")

        results_discordant_regression = rank_regression_test(
            merged_metadata,
            sampling_aggregated_df = median_flux_df,
            group1="PathNonAD-TxAD",
            group2="PathAD-TxNonAD",
            group_col = 'subgroup',
            covariates= ['pmi', 'age_death'],
            categorical_covariates=['batch_x','msex'],
            id_col="projid")
        
        results_discordant_regression["q_value"] = multipletests(
            results_discordant_regression["p_value"],
            method="fdr_bh"
        )[1]
        results_discordant_regression[["compartments_str", "reaction_type", "subsystem", "metabolites"]] = results_discordant_regression["reaction"].apply(get_reaction_info)
        results_discordant_regression = results_discordant_regression.sort_values(by="q_value")
        results_discordant_regression["reaction_name"] = results_discordant_regression["reaction"].map(
            lambda x: model.reactions.get_by_id(x).name
            )
        
        
        output_dir = os.path.join(base_dir, cell_type, "discordant_group")
        os.makedirs(output_dir, exist_ok=True)
        
        results_discordant_regression.to_csv(
            os.path.join(output_dir, "rank_regression_fdr05.csv"),
            index=False
        )
        results_discordant_regression = results_discordant_regression[results_discordant_regression.q_value <= 0.05]
        print(f"Discordant: {len(results_discordant_regression)} reactions at FDR 0.05")

        merged_metadata["path_AD"] = merged_metadata["subgroup"].str.split("-", n=1).str[0] 
        results_PathNon_PathAD_regression = rank_regression_test(
            merged_metadata,
            sampling_aggregated_df = median_flux_df,
            group1="PathNonAD",
            group2="PathAD",
            group_col = 'path_AD',
            covariates= ['pmi', 'age_death'],
            categorical_covariates=['batch_x','msex'],
            id_col="projid")
        
        results_PathNon_PathAD_regression["q_value"] = multipletests(
            results_PathNon_PathAD_regression["p_value"],
            method="fdr_bh"
        )[1]
        
        results_PathNon_PathAD_regression = results_PathNon_PathAD_regression.sort_values(by ='q_value')
        results_PathNon_PathAD_regression["reaction_name"] = results_PathNon_PathAD_regression["reaction"].map(
            lambda x: model.reactions.get_by_id(x).name
            )
        results_PathNon_PathAD_regression[["compartments_str", "reaction_type", "subsystem", "metabolites"]] = results_PathNon_PathAD_regression["reaction"].apply(get_reaction_info)

        output_dir = os.path.join(base_dir, cell_type, "path_ad_path_non_ad_group")
        os.makedirs(output_dir, exist_ok=True)
        
        results_PathNon_PathAD_regression.to_csv(
            os.path.join(output_dir, "rank_regression_fdr05.csv"),
            index=False
        )
        
        results_PathNon_PathAD_regression = results_PathNon_PathAD_regression[results_PathNon_PathAD_regression.q_value <= 0.05]
        print(f"Path: {len(results_PathNon_PathAD_regression)} reactions at FDR 0.05")
        cell_done.append(cell_type)
    except Exception:
        print(f"Error processing cell type: {cell_type}")
        cell_error.append(cell_type)

