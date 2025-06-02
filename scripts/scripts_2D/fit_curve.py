import sys
import os

import numpy as np
import pandas as pd

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parent_dir)
import utils

df = pd.read_csv(snakemake.input[0], index_col=0)
doses = df["dose"].to_numpy()
times = df["time"].to_numpy()
norm_cell_counts = df["norm_cell_count"].to_numpy()
metric = snakemake.config["metric"]
start_time = int(snakemake.config["start_time"])
time = int(snakemake.wildcards["time"])
max_inhibition_time = int(max(times)) # maximal inhibition time in cell count/viability data (minus start_time)
input_parts = snakemake.input[0].split("/")[-1].split("_")
cell_line = input_parts[0]
drug = input_parts[1].split(".")[0]
df_fit_curve_cols = [
    "cell_line", "drug",
    "min_dose", "max_dose",
    "time",
    "alpha", "beta", "gamma", "delta",
    metric
]

dose_differences = np.diff(sorted(set(doses)))
mean_dose_difference = np.mean(dose_differences) # because differences between two consecutive doses are roughly the same in log space
min_dose_bound = np.min(doses) - mean_dose_difference
max_dose_bound = np.max(doses) + mean_dose_difference

param_guesses = [
    [np.max(norm_cell_counts)], # alpha
    [0.001, 0.1], # beta
    [np.min(doses), np.median(doses), np.max(doses)], # gamma
    [np.min(norm_cell_counts)] # delta
]

bounds = (np.array([
    [0, np.inf], # alpha
    [0, np.inf], # beta
    [min_dose_bound, max_dose_bound], # gamma
    [0, np.max(norm_cell_counts)]  # delta
]).T)

time_corrected = time
if time - start_time == max_inhibition_time:
    time_corrected = time - start_time
df_time = df[df["time"] == time_corrected] # train only on datapoints equal to time (minus start_time if time is max_inhibition_time)
doses_time = df_time["dose"].to_numpy()
norm_cell_counts_time = df_time["norm_cell_count"].to_numpy()

best_params, score = utils.fit_model(
    model_function=utils.dose_response_model,
    residual_function=utils.residuals_dose_response_model,
    param_guesses=param_guesses,
    function_input=doses_time,
    function_output=norm_cell_counts_time,
    metric=metric,
    bounds=bounds
)

df_fit_curve = pd.DataFrame([[
    cell_line, drug,
    np.min(doses), np.max(doses),
    time,
    *best_params,
    score
]], columns=df_fit_curve_cols)

df_fit_curve.to_csv(snakemake.output[0])