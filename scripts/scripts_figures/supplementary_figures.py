import sys
import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parent_dir)
import utils

if not os.path.exists("../../visualization_figures"):
    os.makedirs("../../visualization_figures")

colors = px.colors.qualitative.Plotly

data = pd.read_csv("../../output/start_time6h/separate_files_3D/LN229_Drug00.csv", index_col=0) # Only for getting the doses
doses = data["dose"].unique()
xticks = 10**doses
xticks = [f"{x:.3f}" for x in xticks]
a, b, c, d = 1, 1, -1.2, 0.3 # Randomly chosen values
params = (a, b, c, d)
x = np.linspace(-2.67, 1.2, 101)
y = utils.dose_response_model(params, x)
x2 = np.linspace(-2.37, 1, 101)
y2 = utils.dose_response_model(params, x2)

# AUC normalization
fig = go.Figure()
fig.update_layout(
    xaxis_title="Dose (" + u"\u03bc" + "M)",
    yaxis_title="Normalized cell count",
    xaxis = dict(
        tickmode = "array",
        tickvals = xticks,
        ticktext = xticks
    ),
    yaxis1=dict(range=[0, 1.05]),
)
fig.update_xaxes(type="log")

fig.add_trace(go.Scatter(
    x=10**x2,
    y=y2,
    name="Area under the curve",
    line=dict(color="rgba(0, 0, 0, 0)"),
    mode="lines",
    fill="tozeroy",
    fillcolor="rgba(128, 128, 128, 0.5)"
))

fig.add_trace(go.Scatter(
    x=10**x,
    y=y,
    name="Dose-response curve", 
    line=dict(color="black"), 
    mode="lines"
))

fig.add_trace(go.Scatter(
    x=[10**np.min(doses)]*2,
    y=[0, 1.05],
    name="Minimum dose", 
    line=dict(color=colors[0]), 
    mode="lines"
))

fig.add_trace(go.Scatter(
    x=[10**np.max(doses)]*2,
    y=[0, 1.05],
    name="Maximum dose", 
    line=dict(color=colors[1]), 
    mode="lines"
))

fig.add_trace(go.Scatter(
    x=[0, 10**np.max(x)],
    y=[y2[0], y2[0]],
    name="Curve value at minimum dose", 
    line=dict(color=colors[2]), 
    mode="lines"
))

fig.write_image("../../visualization_figures/auc_norm_visualization.pdf", scale=3)

# 4-parameter logistic curve
fig = go.Figure()
fig.update_layout(
    xaxis_title="Dose (" + u"\u03bc" + "M)",
    yaxis_title="Normalized cell count",
    xaxis = dict(
        tickmode = "array",
        tickvals = xticks,
        ticktext = xticks
    ),
    yaxis1=dict(range=[0, 1.05]),
)
fig.update_xaxes(type="log")

fig.add_trace(go.Scatter(
    x=[0, 10**np.max(x)],
    y=[d, d],
    name=u"\u03b4", # delta
    line=dict(color=colors[6]), 
    mode="lines"
))

fig.add_trace(go.Scatter(
    x=[10**c]*2,
    y=[0, 1.05],
    name=u"\u03b3", # gamma
    line=dict(color=colors[5]), 
    mode="lines"
))

fig.add_trace(go.Scatter(
    x=[10**(c - 0.075), 10**(c + 0.125)],
    y=utils.dose_response_model(params, np.array([c - 0.075, c + 0.125])),
    name=u"\u03b2", # beta
    line=dict(color=colors[4]), 
    mode="lines"
))

fig.add_trace(go.Scatter(
    x=[10**(c - 0.075), 10**(c + 0.125)],
    y=utils.dose_response_model(params, np.array([c + 0.125, c + 0.125])),
    name=u"\u03b2", # beta
    line=dict(color=colors[4]), 
    mode="lines",
    fill="tonexty",
    showlegend=False
))

fig.add_trace(go.Scatter(
    x=[10**(c - 0.075)]*2,
    y=utils.dose_response_model(params, np.array([c - 0.075, c + 0.125])),
    name=u"\u03b2", # beta
    line=dict(color=colors[4]), 
    mode="lines",
    showlegend=False
))

fig.add_trace(go.Scatter(
    x=[0, 10**np.max(x)],
    y=[1, 1],
    name=u"\u03b1", # alpha
    line=dict(color=colors[3]), 
    mode="lines"
))

fig.add_trace(go.Scatter(
    x=10**x,
    y=y,
    name="Dose-response curve", 
    line=dict(color="black"), 
    mode="lines"
))

fig.write_image("../../visualization_figures/4pl_visualization.pdf", scale=3)