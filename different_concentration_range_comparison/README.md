Running the Snakemake workflow outputs the file `fit_surface_results_sorted_by_relative_vus_decrease_after_72_daily_False.csv`.

The file `concordance.py` needs the following files, which were obtained as follows:
- `fit_surface_results_sorted_by_relative_vus_decrease_after_72_daily_False_no_min_conc.csv` (the results when omitting the minimum concentration): Run the workflow normally, but in `scripts/preprocessing.py` uncomment lines 85 and 86 (`min_conc = np.nanmin(df_merged["dose"])` and `df_merged = df_merged[df_merged["dose"] != min_conc]`) and add `_no_min_conc` to the outputted file.
- `fit_surface_results_sorted_by_relative_vus_decrease_after_72_daily_False_no_max_conc.csv` (the results when omitting the maximum concentration): Run the workflow normally, but in `scripts/preprocessing.py` uncomment lines 87 and 88 (`max_conc = np.nanmax(df_merged["dose"])` and `df_merged = df_merged[df_merged["dose"] != max_conc]`) and add `_no_max_conc` to the outputted file.
- `fit_surface_results_sorted_by_relative_vus_decrease_after_72_daily_False_no_minmax_conc.csv` (the results when omitting the minimum and maximum concentrations): Run the workflow normally, but in `scripts/preprocessing.py` uncomment lines 85 to 88 (`min_conc = np.nanmin(df_merged["dose"])`, `df_merged = df_merged[df_merged["dose"] != min_conc]`, `max_conc = np.nanmax(df_merged["dose"])`, and `df_merged = df_merged[df_merged["dose"] != max_conc]`) and add `_no_minmax_conc` to the outputted file.



Running the Snakemake 2D workflow outputs the file `fit_curve_results.csv`.

The file `concordance2D.py` needs the following files, which were obtained as follows:
- `fit_curve_results_sorted_by_relative_auc_decrease_after_72_daily_False_no_min_conc.csv` (the results when omitting the minimum concentration): Run the workflow normally, but in `scripts/scripts_2D/preprocessing_celltiterglo.py` uncomment lines 34 and 35 (`min_conc = np.nanmin(df_merged["dose"])` and `df_merged = df_merged[df_merged["dose"] != min_conc]`) and add `_no_min_conc` to the outputted file.
- `fit_curve_results_sorted_by_relative_auc_decrease_after_72_daily_False_no_max_conc.csv` (the results when omitting the maximum concentration): Run the workflow normally, but in `scripts/scripts_2D/preprocessing_celltiterglo.py` uncomment lines 36 and 37 (`max_conc = np.nanmax(df_merged["dose"])` and `df_merged = df_merged[df_merged["dose"] != max_conc]`) and add `_no_max_conc` to the outputted file.
- `fit_curve_results_sorted_by_relative_auc_decrease_after_72_daily_False_no_minmax_conc.csv` (the results when omitting the minimum and maximum concentrations): Run the workflow normally, but in `scripts/scripts_2D/preprocessing_celltiterglo.py` uncomment lines 34 to 37 (`min_conc = np.nanmin(df_merged["dose"])`, `df_merged = df_merged[df_merged["dose"] != min_conc]`, `max_conc = np.nanmax(df_merged["dose"])`, and `df_merged = df_merged[df_merged["dose"] != max_conc]`) and add `_no_minmax_conc` to the outputted file.

Then put the files in this folder and run `curve_results_to_decrease_format.py`.