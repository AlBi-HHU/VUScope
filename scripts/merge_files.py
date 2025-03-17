import pandas as pd

dfs_to_merge = []
for df in snakemake.input:
    dfs_to_merge.append(pd.read_csv(df, index_col=0))

df_merged = pd.concat(dfs_to_merge)
df_merged = df_merged.sort_values(by=["cell_line", "drug", "time"]).reset_index(drop=True)
df_merged.to_csv(snakemake.output[0])