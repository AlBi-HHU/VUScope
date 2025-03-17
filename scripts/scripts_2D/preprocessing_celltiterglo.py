import os

import numpy as np
import pandas as pd

output_dir = snakemake.config["output_dir"]
start_time = int(snakemake.config["start_time"])

dfs_to_merge = []
for input_file in snakemake.input:
    df = pd.read_csv(input_file, index_col=0)
    df_clean = df.dropna().loc[:, ["Patient", "Fluid", "Concentration", "Value"]].sort_values(by=["Fluid", "Concentration"])
    dfs_to_merge.append(df_clean)
df_merged = pd.concat(dfs_to_merge)
df_merged = df_merged.rename(columns={
    "Fluid": "drug",
    "Patient": "cell_line",
    "Concentration": "dose",
    "Value": "norm_cell_count"
})
df_merged[["cell_line", "time"]] = df_merged["cell_line"].str.rsplit("-", n=1, expand=True)
df_merged["time"] = df_merged["time"].astype(int)
df_merged = df_merged[["cell_line", "drug", "dose", "time", "norm_cell_count"]]
df_merged = df_merged.sort_values(by=["cell_line", "drug", "dose", "time"])
df_merged = df_merged.reset_index(drop=True)

# Normalize cell count and log10-transform doses
for (cell_line, drug), df_subset in df_merged.groupby(["cell_line", "drug"]):
     for time, df in df_subset.groupby("time"):
        idx = df.index
        df_merged.loc[idx, "norm_cell_count"] = df["norm_cell_count"]/df["norm_cell_count"].iloc[0]
        df_merged.loc[idx, "dose"] = np.log10(df["dose"])

output_path = f"{output_dir}/start_time{start_time}h/separate_files"
if not os.path.exists(output_path):
    os.makedirs(output_path)

all_cell_lines = df_merged["cell_line"].unique()
for cell_line in all_cell_lines:
    sample_data = df_merged[df_merged["cell_line"] == cell_line]
    all_drugs = sample_data["drug"].unique()
    for drug in all_drugs:
        sample = sample_data[sample_data["drug"] == drug]
        sample = sample.reset_index(drop=True)
        sample.to_csv(f"{output_path}/{cell_line}_{drug}.csv")