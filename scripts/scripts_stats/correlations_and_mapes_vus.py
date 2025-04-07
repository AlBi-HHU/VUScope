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

def pcc_plot(x, y, x_title, y_title, normalized=False):
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
                text=("Normalized " if normalized else "") + f"VUS for {x_title}",
                font=dict(size=25)
            ),
            scaleanchor="y",
            scaleratio=1,
            zeroline=False,
            tickfont=dict(size=20)
        ),
        yaxis=dict(
            title=dict(
                text=("Normalized " if normalized else "") + f"VUS for {y_title}",
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

def ccc_plot(x, y, x_title, y_title, normalized=False):
    ccc = utils.concordance_correlation_coefficient(x, y)
    
    window_pad = 20
    if normalized:
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
                text=("Normalized " if normalized else "") + f"VUS for {x_title}",
                font=dict(size=25)
            ),
            scaleanchor="y",
            scaleratio=1,
            zeroline=False,
            tickfont=dict(size=20)
        ),
        yaxis=dict(
            title=dict(
                text=("Normalized " if normalized else "") + f"VUS for {y_title}",
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

df = pd.read_csv(snakemake.input[0], index_col=0)
metric = "mape"
time_x = int(snakemake.wildcards["time_x"])
time_y = int(snakemake.wildcards["time_y"])
daily_x = snakemake.wildcards["daily_x"] == "True"
daily_y = snakemake.wildcards["daily_y"] == "True"

df_x = df[(df["time"] == time_x) & (df["daily"] == daily_x)].copy()
df_y = df[(df["time"] == time_y) & (df["daily"] == daily_y)].copy()

extrapolated = ""
if time_x < time_y:
    extrapolated = "_extrapolated"
x_title = f"{time_x}h" + (f" (daily)" if daily_x else "")
y_title = f"{time_y}h" + (f" (daily)" if daily_y else "")

fig1, pcc = pcc_plot(df_x["vus" + extrapolated], df_y["vus"], f"{x_title} {extrapolated[1:]}", y_title)
fig1.write_image(snakemake.output[0], scale=3)
fig2, pcc_norm = pcc_plot(df_x["vus_norm" + extrapolated], df_y["vus_norm"], f"{x_title} {extrapolated[1:]}", y_title, normalized=True)
fig2.write_image(snakemake.output[1], scale=3)
fig3, ccc = ccc_plot(df_x["vus" + extrapolated], df_y["vus"], f"{x_title} {extrapolated[1:]}", y_title)
fig3.write_image(snakemake.output[2], scale=3)
fig4, ccc_norm = ccc_plot(df_x["vus_norm" + extrapolated], df_y["vus_norm"], f"{x_title} {extrapolated[1:]}", y_title, normalized=True)
fig4.write_image(snakemake.output[3], scale=3)

score = utils.scoring_function(df_x["vus" + extrapolated], df_y["vus"], metric)
score_norm = utils.scoring_function(df_x["vus_norm" + extrapolated], df_y["vus_norm"], metric)

with open(snakemake.output[4], "w") as f:
    f.write(f"PCC of VUS for {x_title} {extrapolated[1:]} and {y_title}: {pcc}\n".replace("  ", " "))
    f.write(f"CCC of VUS for {x_title} {extrapolated[1:]} and {y_title}: {ccc}\n".replace("  ", " "))
    f.write(f"{metric.upper()} of VUS for {x_title} {extrapolated[1:]} and {y_title}: {score}\n".replace("  ", " "))
    f.write(f"PCC of normalized VUS for {x_title} {extrapolated[1:]} and {y_title}: {pcc_norm}\n".replace("  ", " "))
    f.write(f"CCC of normalized VUS for {x_title} {extrapolated[1:]} and {y_title}: {ccc_norm}\n".replace("  ", " "))
    f.write(f"{metric.upper()} of normalized VUS for {x_title} {extrapolated[1:]} and {y_title}: {score_norm}\n".replace("  ", " "))