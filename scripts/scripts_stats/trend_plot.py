import os

import numpy as np
import pandas as pd

import plotly.graph_objects as go
import plotly.express as px

grayscale = snakemake.config["grayscale"]
color = px.colors.qualitative.Plotly[0]
if grayscale:
    color = "black"

df = pd.read_csv(snakemake.input[0], index_col=0)
df["time"] = df["time"]*24
intermediate_inhibition_time = snakemake.config["intermediate_inhibition_time"]
max_inhibition_time = snakemake.config["max_inhibition_time"]
daily = snakemake.wildcards["daily"] == "True"
df_intermediate_inhibition_time = df[(df["time"] == intermediate_inhibition_time) & (df["daily"] == daily)].copy()
df_max_inhibition_time = df[(df["time"] == max_inhibition_time) & (~df["daily"])].copy()

traces = []

# Connect points between left and middle
for (_, row1), (_, row2) in zip(df_max_inhibition_time.iterrows(), df_intermediate_inhibition_time.iterrows()):
    if (row1["cell_line"], row1["drug"]) == (row2["cell_line"], row2["drug"]):
        traces.append(go.Scatter(
            x=[f"GRIVUS<br>for {max_inhibition_time}h", f"GRIVUS<br>for {intermediate_inhibition_time}h" + (" (daily)" if daily else "")],
            y=[row1["vus_norm"], row2["vus_norm"]],
            mode="lines",
            line=dict(color="black", width=1),
            showlegend=False
        ))

# Connect points between middle and right
for (_, row) in df_intermediate_inhibition_time.iterrows():
    traces.append(go.Scatter(
        x=[f"GRIVUS<br>for {intermediate_inhibition_time}h" + (" (daily)" if daily else ""), f"GRIVUS<br>extrapolated<br>until {max_inhibition_time}h" + (" (daily)" if daily else "")],
        y=[row["vus_norm"], row["vus_norm_extrapolated"]],
        mode="lines",
        line=dict(color="black", width=1),
        showlegend=False
    ))

# left
traces.append(go.Scatter(
    x=[f"GRIVUS<br>for {max_inhibition_time}h"]*len(df_max_inhibition_time["vus_norm"]),
    y=df_max_inhibition_time["vus_norm"],
    mode="markers",
    marker=dict(color=color),
))

# middle
traces.append(go.Scatter(
    x=[f"GRIVUS<br>for {intermediate_inhibition_time}h" + (" (daily)" if daily else "")]*len(df_intermediate_inhibition_time["vus_norm"]),
    y=df_intermediate_inhibition_time["vus_norm"],
    mode="markers",
    marker=dict(color=color),
))

# right
traces.append(go.Scatter(
    x=[f"GRIVUS<br>extrapolated<br>until {max_inhibition_time}h" + (" (daily)" if daily else "")]*len(df_intermediate_inhibition_time["vus_norm_extrapolated"]),
    y=df_intermediate_inhibition_time["vus_norm_extrapolated"],
    mode="markers",
    marker=dict(color=color),
))

layout = go.Layout(
    xaxis=dict(
        tickangle=45,
        zeroline=False
    ),
    yaxis=dict(title="GRIVUS"),
    width=400,
    height=700,
    showlegend=False
)

if grayscale:
    layout.update(
        paper_bgcolor="rgba(255, 255, 255, 1)",
        plot_bgcolor="rgba(255, 255, 255, 1)"
    )

fig = go.Figure(data=traces, layout=layout)
fig.write_image(snakemake.output[0], scale=3)