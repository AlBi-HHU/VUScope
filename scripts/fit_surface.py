import numpy as np
import pandas as pd

import utils

df = pd.read_csv(snakemake.input[0], index_col=0)
doses = df["dose"].to_numpy()
times = df["time"].to_numpy()
norm_cell_counts = df["norm_cell_count"].to_numpy()
metric = snakemake.config["metric"]
start_time = int(snakemake.config["start_time"])/24
time = int(snakemake.wildcards["time"])/24
max_inhibition_time = max(times) # maximal inhibition time in Incucyte data (minus start_time)
daily = snakemake.wildcards["daily"] == "True"
input_parts = snakemake.input[0].split("/")[-1].split("_")
cell_line = input_parts[0]
drug = input_parts[1].split(".")[0]
df_fit_surface_cols = [
    "cell_line", "drug",
    "min_dose", "max_dose",
    "time", "daily",
    "k_alpha", "a_alpha", "k_beta", "k_gamma", "k_delta", "a_delta",
    metric,
    f"{metric}_extrapolated_until_max_inhibition_time"
]

#dose_differences = np.diff(sorted(set(doses)))
#mean_dose_difference = np.mean(dose_differences) # Because differences between two consecutive doses are roughly the same in log space
min_dose_bound = np.min(doses) #- mean_dose_difference # uncomment to marginally improve results, but not deemed necessary
max_dose_bound = np.max(doses) #+ mean_dose_difference # uncomment to marginally improve results, but not deemed necessary

param_guesses = [
    [1], # k_alpha # typical doubling time of in vitro cancer cells is about 1 day
    [100], # a_alpha
    [0.1, 1, 10], # k_beta
    [np.min(doses), np.median(doses), np.max(doses)], # k_gamma
    [0.5], # k_delta
    [10], # a_delta
]

bounds = (np.array([
    [0, 2], # k_alpha # assuming that an in vitro cancer cell with growth-stimulating drugs has a doubling time of at most half a day
    [0, np.inf], # a_alpha
    [0, np.inf], # k_beta
    [min_dose_bound, max_dose_bound], # k_gamma
    [0, 2], # k_delta # assuming that an in vitro cancer cell with growth-stimulating drugs has a doubling time of at most half a day
    [0, np.inf], # a_delta
]).T)

if daily:
    time_corrected = time
    if time - start_time == max_inhibition_time:
        time_corrected = time - start_time
    daily_times = list(np.arange(0, int(time_corrected*24) + 1, 24)/24)
    df_times = df[df["time"].isin(daily_times)] # Train only on datapoints at times = [0, 24, 48, ...]
else:
    time_corrected = time
    if time - start_time == max_inhibition_time:
        time_corrected = time - start_time
    df_times = df[df["time"] <= time_corrected] # Train only on datapoints until time (minus start_time if time is max_inhibition_time)
doses_times = df_times["dose"].to_numpy()
times_times = df_times["time"].to_numpy()
norm_cell_counts_times = df_times["norm_cell_count"].to_numpy()

try:
    best_params, score = utils.fit_model(
        model_function=utils.dose_time_response_model,
        residual_function=utils.residuals_dose_time_response_model,
        param_guesses=param_guesses,
        function_input=(doses_times, times_times),
        function_output=norm_cell_counts_times,
        metric=metric,
        bounds=bounds
    )
except:
    param_guesses = [
        [1], # k_alpha # typical doubling time of in vitro cancer cells is about 1 day
        [100], # a_alpha
        [0.1, 1, 10], # k_beta
        [np.min(doses) + 0.01, np.median(doses), np.max(doses) - 0.01], # k_gamma # sometimes, there is a numerical issue with the trf algorithm
        [0.5], # k_delta
        [10], # a_delta
    ]
    best_params, score = utils.fit_model(
        model_function=utils.dose_time_response_model,
        residual_function=utils.residuals_dose_time_response_model,
        param_guesses=param_guesses,
        function_input=(doses_times, times_times),
        function_output=norm_cell_counts_times,
        metric=metric,
        bounds=bounds
    )

# Calculate score not only on datapoints each 24h, but every available datapoint until time
if daily:
    mask = times <= time_corrected
    doses_until_time = doses[mask]
    times_until_time = times[mask]
    norm_cell_counts_until_time = norm_cell_counts[mask]
    surface_until_time = utils.dose_time_response_model(best_params, (doses_until_time, times_until_time))
    score = utils.scoring_function(norm_cell_counts_until_time, surface_until_time, metric)

# Calculate score between extrapolated function and data points until maximal inhibition time (minus start_time)
score_extrapolated = None
if time_corrected < max_inhibition_time:
    surface_extrapolated = utils.dose_time_response_model(best_params, (doses, times))
    score_extrapolated = utils.scoring_function(norm_cell_counts, surface_extrapolated, metric)

df_fit_surface = pd.DataFrame([[
    cell_line, drug,
    np.min(doses), np.max(doses),
    time, daily,
    *best_params,
    score,
    score_extrapolated
]], columns=df_fit_surface_cols)

df_fit_surface.to_csv(snakemake.output[0])
