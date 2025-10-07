import pandas as pd

for file in ["fit_curve_results", "fit_curve_results_no_min_conc", "fit_curve_results_no_max_conc", "fit_curve_results_no_minmax_conc"]:
    df = pd.read_csv(file + ".csv")
    ending = ""
    if file != "fit_curve_results":
        ending = "_" + "_".join(file.split("_")[3:])

    pivoted = df.pivot_table(
        index=["cell_line", "drug"],
        columns="time",
        values="auc_norm"
    ).reset_index()

    pivoted.columns.name = None
    pivoted = pivoted.rename(columns={
        72: "auc_norm",
        120: "auc_norm_extrapolated"
    })

    pivoted["relative_auc_decrease"] = (pivoted["auc_norm"] - pivoted["auc_norm_extrapolated"])/pivoted["auc_norm"]
    pivoted = pivoted.sort_values(by="relative_auc_decrease", ascending=False).reset_index(drop=True)
    pivoted.to_csv("fit_curve_results_sorted_by_relative_auc_decrease_after_72_daily_False" + ending + ".csv")