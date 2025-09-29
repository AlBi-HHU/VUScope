import numpy as np
import pandas as pd

from scipy.stats import spearmanr

def merge_files(df1, df2):
    df1["pair"] = df1["cell_line"] + "_" + df1["drug"]
    df2["pair"] = df2["cell_line"] + "_" + df2["drug"]

    df1["rank1"] = df1.groupby("cell_line")
    df2["rank2"] = df2.groupby("cell_line")
    
    merged = pd.merge(
        df1[["cell_line", "pair", "rank1"]],
        df2[["cell_line", "pair", "rank2"]],
        on=["cell_line", "pair"]
    )
    
    results = {}
    for cell_line, group in merged.groupby("cell_line"):
        scc, pval = spearmanr(group["rank1"], group["rank2"])
        results[cell_line] = (scc, pval)
        
    return results

relative_vus_decrease = pd.read_csv("fit_surface_results_sorted_by_relative_vus_decrease_after_72_daily_False.csv", index_col=0)
relative_vus_decrease_no_min_conc = pd.read_csv("fit_surface_results_sorted_by_relative_vus_decrease_after_72_daily_False_no_min_conc.csv", index_col=0)
relative_vus_decrease_no_max_conc = pd.read_csv("fit_surface_results_sorted_by_relative_vus_decrease_after_72_daily_False_no_max_conc.csv", index_col=0)
relative_vus_decrease_no_minmax_conc = pd.read_csv("fit_surface_results_sorted_by_relative_vus_decrease_after_72_daily_False_no_minmax_conc.csv", index_col=0)

results_vs_no_min_conc = merge_files(relative_vus_decrease, relative_vus_decrease_no_min_conc)
mean_results_vs_no_min_conc = np.mean([v[0] for v in results_vs_no_min_conc.values()])
print("scc_relative_vs_no_min_conc", results_vs_no_min_conc, "\nmean:", mean_results_vs_no_min_conc)

results_vs_no_max_conc = merge_files(relative_vus_decrease, relative_vus_decrease_no_max_conc)
mean_results_vs_no_max_conc = np.mean([v[0] for v in results_vs_no_max_conc.values()])
print("scc_relative_vs_no_max_conc", results_vs_no_max_conc, "\nmean:", mean_results_vs_no_max_conc)

results_vs_no_minmax_conc = merge_files(relative_vus_decrease, relative_vus_decrease_no_minmax_conc)
mean_results_vs_no_minmax_conc = np.mean([v[0] for v in results_vs_no_minmax_conc.values()])
print("scc_relative_vs_no_minmax_conc", results_vs_no_minmax_conc, "\nmean:", mean_results_vs_no_minmax_conc)

results_no_min_conc_vs_no_max_conc = merge_files(relative_vus_decrease_no_min_conc, relative_vus_decrease_no_max_conc)
mean_results_no_min_conc_vs_no_max_conc = np.mean([v[0] for v in results_no_min_conc_vs_no_max_conc.values()])
print("scc_relative_no_min_conc_vs_no_max_conc", results_no_min_conc_vs_no_max_conc, "\nmean:", mean_results_no_min_conc_vs_no_max_conc)
