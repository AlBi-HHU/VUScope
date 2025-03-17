import os

import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

grayscale = snakemake.config["grayscale"]
color = px.colors.qualitative.Plotly[0]
if grayscale:
    color = "black"

# Distribution of all normalized cell counts to provide the context for metric values
output_dir = snakemake.config["output_dir"]
start_time_incucyte = snakemake.config["start_time_incucyte"]
path = f"{output_dir}/start_time{start_time_incucyte}h/separate_files/"

separate_files = sorted(set(os.listdir(path)) - set([".snakemake_timestamp"]))
dfs_to_merge = []
for df in separate_files:
    dfs_to_merge.append(pd.read_csv(path + df, index_col=0))
df_merged = pd.concat(dfs_to_merge)

norm_cell_counts = list(df_merged["norm_cell_count"])

fig = go.Figure()

fig.add_trace(go.Histogram(
    x=norm_cell_counts,
    xbins=dict(
        start=min(norm_cell_counts),
        end=max(norm_cell_counts),
        size=0.1
    ),
    marker=dict(color=color)
))

fig.update_layout(
    title=f"Distribution of all normalized cell counts (Incucyte)",
    xaxis_title="Normalized cell count",
    yaxis_title="Frequency",
    bargap=0.05
)

if grayscale:
    fig.update_layout(
        paper_bgcolor="rgba(255, 255, 255, 1)",
        plot_bgcolor="rgba(255, 255, 255, 1)"
    )


fig.write_image(snakemake.output[0], scale=3)