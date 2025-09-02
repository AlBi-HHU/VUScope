import numpy as np
import pandas as pd

import utils

def calculate_num_vus(best_params, max_time, min_dose, max_dose, steps=100):
    """
    Calculates the volume under the surface numerically.
    """
    dose_steps = np.linspace(min_dose, max_dose, steps + 1)
    time_steps = np.linspace(0, max_time, steps + 1)
    dose_step = (max_dose - min_dose)/steps
    time_step = max_time/steps

    num_vus_above = 0
    num_vus_below = 0
    norm_factor = 0
    for t in time_steps[:-1]:
        ts = np.array([t, t + time_step, t, t + time_step])
        for d in dose_steps[:-1]:
            ds = np.array([d, d, d + dose_step, d + dose_step])
            approx = utils.dose_time_response_model_grivus(best_params, (ds, ts)) # remove "_grivus" for normal vus
            approx_above = np.max(approx)
            approx_below = np.min(approx)
            num_vus_above += approx_above*(time_step*dose_step)
            num_vus_below += approx_below*(time_step*dose_step)
        approx_at_min_dose = utils.dose_time_response_model_grivus(best_params, (np.array([min_dose, min_dose]), np.array([t, t + time_step]))) # remove "_grivus" for normal vus
        approx_at_min_dose_above = np.max(approx_at_min_dose)
        approx_at_min_dose_below = np.min(approx_at_min_dose)
        num_auc_at_min_dose_above = approx_at_min_dose_above*time_step
        num_auc_at_min_dose_below = approx_at_min_dose_below*time_step
        norm_factor += (num_auc_at_min_dose_above + num_auc_at_min_dose_below)/2
    norm_factor *= (max_dose - min_dose)

    vus = (num_vus_above + num_vus_below)/2
    vus_norm = vus/norm_factor
    
    return vus_norm

dfs_to_merge = []
for df in snakemake.input:
    dfs_to_merge.append(pd.read_csv(df, index_col=0))
df_merged = pd.concat(dfs_to_merge).reset_index(drop=True)
times = sorted(df_merged["time"].unique())
dailies = df_merged["daily"].unique()
start_time = snakemake.config["start_time"]/24
max_inhibition_time = max(times) - start_time # Maximal inhibition time in IncuCyte data (minus start_time)
extrapolate_until_time = snakemake.config["extrapolate_until_time"]/24

dfs_surface_fit = []
for i in times:
    for j in dailies:
        df_surface_fit = df_merged[(df_merged["time"] == i) & (df_merged["daily"] == j)].copy()
        best_params = np.array(df_surface_fit.iloc[0, 6:-2])
        
        # Calculate VUS until time
        # If max_time is greater than max_inhibition_time, we calculate the VUS until max_inhibition_time
        min_dose, max_dose = df_surface_fit.iloc[0, 2:4]
        time = df_surface_fit.iloc[0, 4]
        if time > max_inhibition_time:
            time = max_inhibition_time
        df_surface_fit["vus_norm"] = calculate_num_vus(best_params, time, min_dose, max_dose)
        
        # Calculate VUS of extrapolated surface until the value of extrapolate_until_time (minus start_time) in config.yaml
        max_time = extrapolate_until_time - start_time
        df_surface_fit["vus_norm_extrapolated"] = calculate_num_vus(best_params, max_time, min_dose, max_dose)

        dfs_surface_fit.append(df_surface_fit)

df_surface_fits = pd.concat(dfs_surface_fit)
df_surface_fits.to_csv(snakemake.output[0])