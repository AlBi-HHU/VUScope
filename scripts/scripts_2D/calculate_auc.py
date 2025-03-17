import sys
import os

import numpy as np
import pandas as pd

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parent_dir)
import utils

dfs_to_merge = []
for df in snakemake.input:
    dfs_to_merge.append(pd.read_csv(df, index_col=0))
df_merged = pd.concat(dfs_to_merge).reset_index(drop=True)
times = sorted(df_merged["time"].unique())
start_time = snakemake.config["start_time"]

dfs_curve_fit = []
for i in times:
    df_curve_fit = df_merged[df_merged["time"] == i].copy()
    best_params = np.array(df_curve_fit.iloc[0, 5:-1])
    
    # Calculate AUC at time
    min_dose, max_dose = df_curve_fit.iloc[0, 2:4]
    df_curve_fit["auc"], df_curve_fit["auc_norm"] = utils.dose_response_model_auc(best_params, min_dose, max_dose)
    dfs_curve_fit.append(df_curve_fit)

df_curve_fits = pd.concat(dfs_curve_fit)
df_curve_fits.to_csv(snakemake.output[0])