import sys
import os

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, linregress
import plotly.graph_objects as go
import plotly.express as px

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parent_dir)
import utils

grayscale = snakemake.config["grayscale"]
color = px.colors.qualitative.Plotly[0]
if grayscale:
    color = "black"

def pcc_plot(x, y, x_title, y_title):
    pcc, p = pearsonr(x, y)
    
    slope, intercept, _, _, _ = linregress(x, y)
    line_y = slope*x + intercept

    scatter = go.Scatter(
        x=x,
        y=y,
        mode="markers",
        marker=dict(color=color)
    )
    line = go.Scatter(
        x=x,
        y=line_y,
        mode="lines",
        line=dict(color=color)
    )
    layout = go.Layout(
        title=dict(
            text=f"Pearson Correlation Coefficient<br>(r = {pcc:.3f}, p-value = {p:.1e})",
            font=dict(size=30)
        ),
        xaxis=dict(
            title=dict(
                text=f"Normalized VUS for {x_title}",
                font=dict(size=25)
            ),
            scaleanchor="y",
            scaleratio=1,
            zeroline=False,
            tickfont=dict(size=20)
        ),
        yaxis=dict(
            title=dict(
                text=f"Normalized AUC for {y_title}",
                font=dict(size=25)
            ),
            tickfont=dict(size=20)
        ),
        width=800,
        height=800,
        showlegend=False
    )
    if grayscale:
        layout.update(
            paper_bgcolor="rgba(255, 255, 255, 1)",
            plot_bgcolor="rgba(255, 255, 255, 1)"
        )
    fig = go.Figure(data=[scatter, line], layout=layout)
    return fig, pcc

def ccc_plot(x, y, x_title, y_title):
    ccc = utils.concordance_correlation_coefficient(x, y)
    window_pad = 0.02
    minimum = min(x.min(), y.min()) - window_pad
    maximum = max(x.max(), y.max()) + window_pad
    
    scatter = go.Scatter(
        x=x,
        y=y,
        mode="markers",
        marker=dict(color=color)
    )
    line = go.Scatter(
        x=[minimum, maximum],
        y=[minimum, maximum],
        mode="lines",
        line=dict(color=color)
    )
    layout = go.Layout(
        title=dict(
            text=f"Concordance Correlation Coefficient<br>(CCC = {ccc:.3f})",
            font=dict(size=30)
        ),
        xaxis=dict(
            title=dict(
                text=f"Normalized VUS for {x_title}",
                font=dict(size=25)
            ),
            scaleanchor="y",
            scaleratio=1,
            zeroline=False,
            tickfont=dict(size=20)
        ),
        yaxis=dict(
            title=dict(
                text=f"Normalized AUC for {y_title}",
                font=dict(size=25)
            ),
            tickfont=dict(size=20)
        ),
        width=800,
        height=800,
        showlegend=False
    )
    if grayscale:
        layout.update(
            paper_bgcolor="rgba(255, 255, 255, 1)",
            plot_bgcolor="rgba(255, 255, 255, 1)"
        )
    fig = go.Figure(data=[scatter, line], layout=layout)
    return fig, ccc

df_incucyte = pd.read_csv(snakemake.input[0], index_col=0)
df_incucyte["time"] = df_incucyte["time"]*24
df_celltiterglo = pd.read_csv(snakemake.input[1], index_col=0)
max_inhibition_time = snakemake.config["max_inhibition_time"]
metric = "mape"
time_x = int(snakemake.wildcards["time_x"])
daily_x = snakemake.wildcards["daily_x"] == "True"
time_y = int(snakemake.wildcards["time_y"])

df_x = df_incucyte[(df_incucyte["time"] == time_x) & (df_incucyte["daily"] == daily_x)].copy()
df_y = df_celltiterglo[df_celltiterglo["time"] == time_y].copy()

extrapolated = ""
if time_x < time_y:
    extrapolated = "_extrapolated"
x_title = f"{time_x}h Incucyte" + (f" (daily)" if daily_x else "")
y_title = f"{time_y}h CellTiter-Glo"

fig, pcc = pcc_plot(df_x["vus_norm" + extrapolated], df_y["auc_norm"], f"{x_title} {extrapolated[1:]}", y_title)
fig.write_image(snakemake.output[0], scale=3)
fig, ccc = ccc_plot(df_x["vus_norm" + extrapolated], df_y["auc_norm"], f"{x_title} {extrapolated[1:]}", y_title)
fig.write_image(snakemake.output[1], scale=3)

score = utils.scoring_function(df_x["vus_norm" + extrapolated], df_y["auc_norm"], metric)

with open(snakemake.output[2], "w") as f:
    f.write(f"PCC of normalized VUS for {x_title} {extrapolated[1:]} and AUC for {y_title}: {pcc}\n".replace("  ", " "))
    f.write(f"CCC of normalized VUS for {x_title} {extrapolated[1:]} and AUC for {y_title}: {ccc}\n".replace("  ", " "))
    f.write(f"{metric.upper()} of normalized VUS for {x_title} {extrapolated[1:]} and AUC for {y_title}: {score}\n".replace("  ", " "))