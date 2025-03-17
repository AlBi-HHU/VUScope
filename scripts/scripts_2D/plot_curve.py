import sys
import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parent_dir)
import utils

steps_curve = 20
colors = ["grey", "black"]
separate_file_sample = pd.read_csv(snakemake.input[0], index_col=0)
dfs_to_merge = []
for df in snakemake.input[1:]:
    dfs_to_merge.append(pd.read_csv(df, index_col=0))
df_merged = pd.concat(dfs_to_merge).reset_index(drop=True)
start_time = snakemake.config["start_time"]
    
fig = go.Figure()
doses = separate_file_sample["dose"]
xticks = 10**doses
xticks = [f"{x:.3f}" for x in xticks]

for i in range(df_merged.shape[0]):
    cell_line, drug = df_merged.iloc[i, :2]
    min_dose, max_dose = df_merged.iloc[i, 2:4]
    time = df_merged.iloc[i, 4]
    best_params = np.array(df_merged.iloc[i, 5:-1])
    
    # Draw data points
    time_corrected = time
    if time - start_time == max(separate_file_sample["time"]):
        time_corrected = time - start_time
    mask = separate_file_sample["time"] == time_corrected
    norm_cell_count_norm = separate_file_sample["norm_cell_count"][mask]/(np.array(separate_file_sample["norm_cell_count"][mask])[0]) # Normalize to have first data point at y=0, which is common for 2D curve fits
    fig.add_trace(go.Scatter(
        x=separate_file_sample["dose"][mask],
        y=norm_cell_count_norm,
        mode="markers",
        marker=dict(color=colors[i]),
        opacity=0.7 if i == 0 else 1,
        showlegend=False
    ))
    
    # Draw curve    
    dose_steps = np.linspace(min_dose, max_dose, len(separate_file_sample["dose"].unique())*steps_curve + 1)
    norm_cell_count_steps = np.array(utils.dose_response_model(best_params, dose_steps))
    norm_cell_count_steps_norm = norm_cell_count_steps/norm_cell_count_steps[0]
    fig.add_trace(go.Scatter(
        x=dose_steps,
        y=norm_cell_count_steps_norm,
        name=f"{time}h (" + u"\u03b3" + f" = {best_params[2]:.3f})",
        mode="lines",
        line=dict(color=colors[i]),
        opacity=0.7 if i == 0 else 1
    ))

fig.update_layout(
    title=dict(
        text=f"Cell line {cell_line} with drug {drug}",
        font=dict(size=30)
    ),
    xaxis=dict(
        title=dict(
            text="Dose (" + u"\u03bc" + "M)",
            font=dict(size=25)
        ),
        tickmode="array",
        tickvals=doses,
        ticktext=xticks,
        zeroline=False,
        tickfont=dict(size=20)
    ),
    yaxis=dict(
        title=dict(
            text="Normalized cell count",
            font=dict(size=25)
        ),
        tickfont=dict(size=20)
    ),
    legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1,
        xanchor="center",
        x=0.5,
        font=dict(size=25)
    ),
    paper_bgcolor="rgba(255, 255, 255, 1)",
    plot_bgcolor="rgba(255, 255, 255, 1)",
    width=800,
    height=800
)

fig.write_image(snakemake.output[0], scale=3)