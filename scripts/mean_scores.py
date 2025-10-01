import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

import utils

def mean_score(df_merged, times, metric, daily=None, max_inhibition_time=None):
    if daily is None:
        mean_score = df_merged[df_merged["time"] == time][metric].mean()
        std = df_merged[df_merged["time"] == time][metric].std()
        score_string = f"Mean {metric.upper()} from curve to datapoints at {time}h: " + str(mean_score) + " " + u"\u00B1" + " " + str(std) + "\n"
    else:
        mean_score = df_merged[(df_merged["time"] == time) & (df_merged["daily"] == d)][metric].mean()
        std = df_merged[(df_merged["time"] == time) & (df_merged["daily"] == d)][metric].std()
        score_string = f"Mean {metric.upper()} from surface to datapoints for {time}h" + (" (daily)" if d else "") + ": " + str(mean_score) + " " + u"\u00B1" + " " + str(std) + "\n"
        if time < max_inhibition_time:
            mean_score_extrapolated_until_max_inhibition_time = df_merged[df_merged["daily"] == d][f"{metric}_extrapolated_until_max_inhibition_time"].mean()
            std_extrapolated_until_max_inhibition_time = df_merged[df_merged["daily"] == d][f"{metric}_extrapolated_until_max_inhibition_time"].std()
            score_string += f"Mean {metric.upper()} from surface to datapoints for {time}h" + (" (daily)" if d else "") + f" extrapolated until {max_inhibition_time}h: " + str(mean_score_extrapolated_until_max_inhibition_time) + " " + u"\u00B1" + " " + str(std_extrapolated_until_max_inhibition_time) + "\n"
    return score_string

df_merged = pd.read_csv(snakemake.input[0], index_col=0)
output_dir = snakemake.config["output_dir"]
metric = snakemake.config["metric"]
times = np.array(snakemake.config["times"])/24

if "/2D/" in snakemake.input[0]:
    with open(snakemake.output[0], "w") as f:
        for time in times:
            score_string = mean_score(df_merged, time, metric)
            f.write(score_string)
else:
    max_inhibition_time = max(times) # maximal inhibition time in Incucyte data
    daily = snakemake.config["daily"]

    with open(snakemake.output[0], "w") as f:
        for d in daily:
            for time in times:
                score_string = mean_score(df_merged, time, metric, d, max_inhibition_time)
                f.write(score_string)