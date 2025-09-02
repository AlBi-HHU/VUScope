import pandas as pd

df_merged = pd.read_csv(snakemake.input[0], index_col=0)
time = int(snakemake.wildcards["time"])/24
daily = snakemake.wildcards["daily"] == "True"
df_surface_fit = df_merged[(df_merged["time"] == time) & (df_merged["daily"] == daily)].copy()

df_final = df_surface_fit.loc[:, ["cell_line", "drug", "vus_norm", "vus_norm_extrapolated"]]
df_final["absolute_vus_decrease"] = df_final["vus_norm"] - df_final["vus_norm_extrapolated"]
df_final["relative_vus_decrease"] = (df_final["vus_norm"] - df_final["vus_norm_extrapolated"])/df_final["vus_norm"]
if snakemake.wildcards["sort_by_vus_decrease"] == "relative":
    sort_by_vus_decrease = "relative_vus_decrease"
elif snakemake.wildcards["sort_by_vus_decrease"] == "absolute":
    sort_by_vus_decrease = "absolute_vus_decrease"
else:
    print('Invalid sorting option for VUS decrease in config.yaml, defaulted to "relative"')
    sort_by_vus_decrease = "relative_vus_decrease"
df_final = df_final.sort_values(by=sort_by_vus_decrease, ascending=False).reset_index(drop=True)
df_final.to_csv(snakemake.output[0])