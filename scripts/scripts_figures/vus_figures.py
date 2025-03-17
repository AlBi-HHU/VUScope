import sys
import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parent_dir)
import utils

if not os.path.exists("../../visualization_figures"):
    os.makedirs("../../visualization_figures")

steps_surface = 5
colors = px.colors.qualitative.Plotly
colorscale = "Plotly3"
custom_colorscale = [
    [0.0, "grey"],
    [1.0, "grey"]
]

separate_file_sample = pd.read_csv("../../output/start_time6h/separate_files/UW228-3_Drug11.csv", index_col=0)
df_merged = pd.read_csv("../../output/start_time6h/3D/extrapolate_until_120h/tmp/fit_surface/UW228-3_Drug11_time72_dailyTrue.csv", index_col=0)
df_merged2 = pd.read_csv("../../output/start_time6h/3D/extrapolate_until_120h/tmp/fit_surface/UW228-3_Drug11_time120_dailyTrue.csv", index_col=0)
doses = separate_file_sample["dose"].unique()
xticks = 10**doses
xticks = [f"{x:.3f}" for x in xticks]

cell_line, drug = df_merged.iloc[0, :2]
min_dose, max_dose = df_merged.iloc[0, 2:4]
max_time = separate_file_sample["time"].max()
best_params = np.array(df_merged.iloc[0, 6:-2])
dose_steps = np.linspace(min_dose, max_dose, len(separate_file_sample["dose"].unique())*steps_surface + 1)
time_steps = np.linspace(0, max_time, len(separate_file_sample["time"].unique())*steps_surface + 1)
norm_cell_count_steps = np.array([utils.dose_time_response_model(best_params, (dose_steps, t)) for t in time_steps])

z_72h = utils.dose_time_response_model(best_params, (dose_steps, 72))
z_120h = utils.dose_time_response_model(best_params, (dose_steps, max_time))
value_at_min_dose_and_72h = utils.dose_time_response_model(best_params, (min_dose, 72))
value_at_min_dose_and_120h = utils.dose_time_response_model(best_params, (min_dose, max_time))

# Determine maximum height (z-axis value) to use for all plots
max_z = 0

best_params = np.array(df_merged.iloc[0, 6:-2])
norm_cell_count_steps = np.array([utils.dose_time_response_model(best_params, (dose_steps, time_step)) for time_step in time_steps])
max_z = np.max([max_z, np.max(norm_cell_count_steps)])

best_params2 = np.array(df_merged2.iloc[0, 6:-2])
norm_cell_count_steps2 = np.array([utils.dose_time_response_model(best_params2, (dose_steps, time_step)) for time_step in time_steps])
max_z = np.max([max_z, np.max(norm_cell_count_steps2)])

# VUS
fig = go.Figure()

fig.update_layout(
    scene=dict(
        #xaxis_title="Dose (" + u"\u03bc" + "M)",
        #yaxis_title="Time (h)",
        #zaxis_title="Normalized cell count",
        xaxis_title="",
        yaxis_title="",
        zaxis_title="",
        xaxis=dict(
            tickmode="array",
            tickvals=doses,
            ticktext=xticks,
            zeroline=False
        ),
        zaxis=dict(
            range=[0, max_z + 0.05],
            tickvals=[""] + list(np.arange((max_z + 0.05) + 1))[1:],
            autorange=False
        )
    ),
    scene_camera={"eye": {"x": 1.8, "y": -1.8, "z": .5}}
)

fig.add_trace(go.Surface(
    x=dose_steps,
    y=time_steps,
    z=norm_cell_count_steps,
    colorscale=colorscale,
    name=f"{cell_line}<br>{drug}",
    opacity=0.5,
    showscale=False
))

fig.add_trace(go.Surface(
    x=dose_steps,
    y=time_steps,
    z=np.zeros_like(norm_cell_count_steps),
    colorscale=custom_colorscale,
    opacity=0.5,
    showscale=False
))

fig.add_trace(go.Surface(
    x=[np.min(dose_steps), np.max(dose_steps)], 
    y=[0, 0],
    z=[[0, 0], [1, 1]],
    colorscale=custom_colorscale,
    opacity=0.5,
    showscale=False
))

fig.add_trace(go.Surface(
    x=[np.min(dose_steps), np.max(dose_steps)],
    y=[0, 0],
    z=[[0, 0], [1, 1]],
    colorscale=custom_colorscale,
    opacity=0.5,
    showscale=False
))

fig.add_trace(go.Surface(
    x=dose_steps,
    y=[max_time, max_time],
    z=[[0]*len(dose_steps), utils.dose_time_response_model(best_params, (dose_steps, max_time))],
    colorscale=custom_colorscale,
    opacity=0.5,
    showscale=False
))

fig.add_trace(go.Surface(
    x=np.array([np.min(dose_steps), np.min(dose_steps)]), 
    y=time_steps,
    z=[[0, utils.dose_time_response_model(best_params, (np.min(dose_steps), t))] for t in time_steps],
    colorscale=custom_colorscale,
    opacity=0.5,
    showscale=False
))

fig.add_trace(go.Surface(
    x=np.array([np.max(dose_steps), np.max(dose_steps)]), 
    y=time_steps,
    z=[[0, utils.dose_time_response_model(best_params, (np.max(dose_steps), t))] for t in time_steps],
    colorscale=custom_colorscale,
    opacity=0.5,
    showscale=False
))
fig.add_trace(go.Scatter3d(
    x=[min_dose, min_dose],
    y=[max_time, max_time],
    z=[0, value_at_min_dose_and_120h],
    mode="lines",
    opacity=0.2,
    line=dict(color="grey")
))

fig.write_image("../../visualization_figures/vus_visualization.pdf", scale=3)

# Numerical VUS
steps = 5
dose_steps_num_vus = np.linspace(min_dose, max_dose, steps + 1)
time_steps_num_vus = np.linspace(0, max_time, steps + 1)
dose_step = (max_dose - min_dose)/steps
time_step = max_time/steps

# Numerical VUS above
fig = go.Figure()

fig.add_trace(go.Surface(
    x=dose_steps,
    y=time_steps,
    z=norm_cell_count_steps,
    colorscale=colorscale,
    opacity=0.7,
    showscale=False
))

for d in dose_steps_num_vus[:-1]:
    ds = np.array([d, d, d + dose_step, d + dose_step])
    for t in time_steps_num_vus[:-1]:
        ts = np.array([t, t + time_step, t, t + time_step])
        approx = utils.dose_time_response_model(best_params, (ds, ts))
        approx_above = np.max(approx)
        
        fig.add_trace(go.Surface(
            x=[d, d + dose_step],
            y=[t, t + time_step],
            z=[[approx_above, approx_above]]*2,
            colorscale=custom_colorscale,
            showscale=False,
            opacity=0.5
        ))
        
        fig.add_trace(go.Surface(
            x=[d, d + dose_step],
            y=[t, t],
            z=[[approx_above, approx_above], [0, 0]],
            colorscale=custom_colorscale,
            showscale=False,
            opacity=0.5
        ))
        
        fig.add_trace(go.Surface(
            x=[d + dose_step, d + dose_step],
            y=[t, t + time_step],
            z=[[0, approx_above], [0, approx_above]],
            colorscale=custom_colorscale,
            showscale=False,
            opacity=0.3
        ))
    
fig.update_layout(
    scene=dict(
        #xaxis_title="Dose (" + u"\u03bc" + "M)",
        #yaxis_title="Time (h)",
        #zaxis_title="Normalized cell count",
        xaxis_title="",
        yaxis_title="",
        zaxis_title="",
        xaxis=dict(
            tickmode="array",
            tickvals=doses,
            ticktext=xticks,
            zeroline=False
        ),
        zaxis=dict(
            range=[0, max_z + 0.05],
            tickvals=[""] + list(np.arange((max_z + 0.05) + 1))[1:],
            autorange=False
        )
    ),
    scene_camera={"eye": {"x": 1.8, "y": -1.8, "z": .5}}
)

fig.write_image("../../visualization_figures/numerical_vus_above.pdf", scale=3)

# Numerical VUS below
fig = go.Figure()

fig.add_trace(go.Surface(
    x=dose_steps,
    y=time_steps,
    z=norm_cell_count_steps,
    colorscale=colorscale,
    opacity=0.7,
    showscale=False
))

for d in dose_steps_num_vus[:-1]:
    ds = np.array([d, d, d + dose_step, d + dose_step])
    for t in time_steps_num_vus[:-1]:
        ts = np.array([t, t + time_step, t, t + time_step])
        approx = utils.dose_time_response_model(best_params, (ds, ts))
        approx_below = np.min(approx)
        
        fig.add_trace(go.Surface(
            x=[d, d + dose_step],
            y=[t, t + time_step],
            z=[[approx_below, approx_below]]*2,
            colorscale=custom_colorscale,
            showscale=False,
            opacity=0.5
        ))
        
        fig.add_trace(go.Surface(
            x=[d, d + dose_step],
            y=[t, t],
            z=[[approx_below, approx_below], [0, 0]],
            colorscale=custom_colorscale,
            showscale=False,
            opacity=0.5
        ))
        
        fig.add_trace(go.Surface(
            x=[d + dose_step, d + dose_step],
            y=[t, t + time_step],
            z=[[0, approx_below], [0, approx_below]],
            colorscale=custom_colorscale,
            showscale=False,
            opacity=0.3
        ))

fig.update_layout(
    scene=dict(
        #xaxis_title="Dose (" + u"\u03bc" + "M)",
        #yaxis_title="Time (h)",
        #zaxis_title="Normalized cell count",
        xaxis_title="",
        yaxis_title="",
        zaxis_title="",
        xaxis=dict(
            tickmode="array",
            tickvals=doses,
            ticktext=xticks,
            zeroline=False
        ),
        zaxis=dict(
            range=[0, max_z + 0.05],
            tickvals=[""] + list(np.arange((max_z + 0.05) + 1))[1:],
            autorange=False
        )
    ),
    scene_camera={"eye": {"x": 1.8, "y": -1.8, "z": .5}}
)

fig.write_image("../../visualization_figures/numerical_vus_below.pdf", scale=3)

# VUS normalization
"""
fig = make_subplots(
    rows=2, cols=3,
    specs=[
        [{}, {"rowspan": 2, "colspan": 2, "type": "scatter3d"}, None],
        [{}, None, None]
    ],
    subplot_titles=("At 72h", None, "At 120h", None)
)

fig["layout"]["xaxis"]["title"] = "Dose (" + u"\u03bc" + "M)"
fig["layout"]["xaxis2"]["title"] = "Dose (" + u"\u03bc" + "M)"
fig["layout"]["yaxis"]["title"] = "Normalized cell count"
fig["layout"]["yaxis2"]["title"] = "Normalized cell count"

fig.update_layout(
    xaxis=dict(
        tickmode="array",
        tickvals=doses,
        ticktext=xticks,
        zeroline=False
    ),
    xaxis2=dict(
        tickmode="array",
        tickvals=doses,
        ticktext=xticks,
        zeroline=False
    ),
    width=1000,
    height=700,
)

# Upper left
fig.add_trace(go.Scatter(
    x=dose_steps,
    y=z_72h,
    name="Dose-response curve",
    line=dict(color="grey", width=3),
    mode="lines",
    showlegend=False
), row=1, col=1)

fig.add_trace(go.Scatter(
    x=[min_dose, min_dose],
    y=[0, value_at_min_dose_and_72h],
    name="Minimum dose",
    mode="lines",
    line=dict(color=colors[0], width=3),
    showlegend=False
), row=1, col=1)

fig.add_trace(go.Scatter(
    x=[max_dose, max_dose],
    y=[0, value_at_min_dose_and_72h],
    name="Maximum dose",
    mode="lines",
    line=dict(color=colors[1], width=3),
    showlegend=False
), row=1, col=1)

fig.add_trace(go.Scatter(
    x=[min_dose, max_dose],
    y=[value_at_min_dose_and_72h, value_at_min_dose_and_72h],
    name="Value at minimum dose at 72h",
    mode="lines",
    line=dict(color=colors[7], width=3),
    showlegend=False
), row=1, col=1)

fig.add_trace(go.Scatter(
    x=[min_dose, max_dose],
    y=[0, 0],
    name="Maximum dose minus minimum dose",
    mode="lines",
    line=dict(color=colors[5], width=3),
    showlegend=False
), row=1, col=1)

# Lower left
fig.add_trace(go.Scatter(
    x=dose_steps,
    y=z_120h,
    name="Dose-response curve",
    line=dict(color="black", width=3),
    mode="lines",
    showlegend=False
), row=2, col=1)

fig.add_trace(go.Scatter(
    x=[min_dose, min_dose],
    y=[0, value_at_min_dose_and_120h],
    name="Minimum dose",
    mode="lines",
    line=dict(color=colors[0]),
    showlegend=False
), row=2, col=1)

fig.add_trace(go.Scatter(
    x=[max_dose, max_dose],
    y=[0, value_at_min_dose_and_120h],
    name="Maximum dose",
    mode="lines",
    line=dict(color=colors[1]),
    showlegend=False
), row=2, col=1)

fig.add_trace(go.Scatter(
    x=[min_dose, max_dose],
    y=[value_at_min_dose_and_120h, value_at_min_dose_and_120h],
    name="Value at minimum dose at 120h",
    mode="lines",
    line=dict(color=colors[2]),
    showlegend=False
), row=2, col=1)

fig.add_trace(go.Scatter(
    x=[min_dose, max_dose],
    y=[0, 0],
    name="Maximum dose minus minimum dose",
    mode="lines",
    line=dict(color=colors[5]),
    showlegend=False
), row=2, col=1)

fig.update_layout(
    yaxis1=dict(range=[0, 9]),
    yaxis2=dict(range=[0, 9])
)

fig.update_layout(
    scene=dict(
        #xaxis_title="Dose (" + u"\u03bc" + "M)",
        #yaxis_title="Time (h)",
        #zaxis_title="Normalized cell count",
        xaxis_title="",
        yaxis_title="",
        zaxis_title="",
        xaxis=dict(
            tickmode="array",
            tickvals=doses,
            ticktext=xticks,
            zeroline=False
        ),
        zaxis=dict(
            range=[0, max_z + 0.05],
            tickvals=[""] + list(np.arange((max_z + 0.05) + 1))[1:],
            autorange=False
        )
    ),
    scene_camera={"eye": {"x": 1.8, "y": -1.8, "z": .5}}
)
"""

fig = go.Figure()
# Right
"""
fig.add_trace(go.Scatter3d(
    x=dose_steps,
    y=np.array([72]*101),
    z=z_72h,
    mode="lines",
    line=dict(color="grey", width=3),
    name="Dose-response curve at 72h"
), row=1, col=2)

fig.add_trace(go.Scatter3d(
    x=dose_steps,
    y=np.array([max_time]*101),
    z=z_120h,
    mode="lines",
    line=dict(color="black", width=3),
    name="Dose-response curve at 120h"
), row=1, col=2)

fig.add_trace(go.Scatter3d(
    x=[min_dose, min_dose],
    y=[max_time, max_time],
    z=[0, value_at_min_dose_and_120h],
    mode="lines",
    line=dict(color=colors[0], width=3),
    name="Minimum dose"
), row=1, col=2)

fig.add_trace(go.Scatter3d(
    x=[min_dose, min_dose],
    y=[72, 72],
    z=[0, value_at_min_dose_and_72h],
    mode="lines",
    line=dict(color=colors[0], width=3),
    showlegend=False
), row=1, col=2)

fig.add_trace(go.Scatter3d(
    x=[max_dose, max_dose],
    y=[72, 72],
    z=[0, value_at_min_dose_and_72h],
    mode="lines",
    line=dict(color=colors[1], width=3),
    name="Maximum dose"
), row=1, col=2)

fig.add_trace(go.Scatter3d(
    x=[max_dose, max_dose],
    y=[max_time, max_time],
    z=[0, value_at_min_dose_and_120h],
    mode="lines",
    line=dict(color=colors[1], width=3),
    showlegend=False
), row=1, col=2)

fig.add_trace(go.Scatter3d(
    x=[min_dose, max_dose],
    y=[72, 72],
    z=[0, 0],
    name="Maximum dose minus minimum dose",
    mode="lines",
    line=dict(color=colors[5], width=3)
), row=1, col=2)

fig.add_trace(go.Scatter3d(
    x=[min_dose, max_dose],
    y=[max_time, max_time],
    z=[0, 0],
    name="Maximum dose minus minimum dose",
    mode="lines",
    line=dict(color=colors[5], width=3),
    showlegend=False
), row=1, col=2)

fig.add_trace(go.Scatter3d(
    x=[min_dose, max_dose],
    y=[72, 72],
    z=[value_at_min_dose_and_72h, value_at_min_dose_and_72h],
    mode="lines",
    line=dict(color=colors[7], width=3),
    name="Value at minimum dose at 72h"
), row=1, col=2)

fig.add_trace(go.Scatter3d(
    x=[min_dose, max_dose],
    y=[max_time, max_time],
    z=[value_at_min_dose_and_120h, value_at_min_dose_and_120h],
    mode="lines",
    line=dict(color=colors[2], width=3),
    name="Value at minimum dose at 120h"
), row=1, col=2)

t = np.linspace(0, max_time, 101)
z2 = utils.dose_time_response_model(best_params, (min_dose, t))
fig.add_trace(go.Scatter3d(
    x=np.array([min_dose]*101),
    y=t,
    z=z2,
    mode="lines",
    line=dict(color=colors[4], width=3),
    name="Time-response curve at min dose"
), row=1, col=2)
"""
fig.add_trace(go.Surface(
    x=dose_steps,
    y=time_steps,
    z=norm_cell_count_steps,
    name=f"{cell_line}<br>{drug}",
    opacity=0.5,
    colorscale=colorscale,
    showscale=False
)) #, row=1, col=2)

dose_steps2 = np.array([min_dose, max_dose])
norm_cell_count_steps2 = np.array([utils.dose_time_response_model(best_params, (np.array([min_dose, min_dose]), t)) for t in time_steps])
fig.add_trace(go.Surface(
    x=dose_steps2,
    y=time_steps,
    z=norm_cell_count_steps2,
    name="Normalization",
    opacity=0.5,
    colorscale="greens",
    showscale=False
)) #, row=1, col=2)

"""
fig.update_layout(
    legend=dict(
        yanchor="top",
        y=1,
        xanchor="right",
        x=0.8
    ),
    scene=dict(
        domain=dict(
            x=[0.2, 1],
            y=[0, 0.9]
        )
    )
)
"""

fig.update_layout(
    scene=dict(
        #xaxis_title="Dose (" + u"\u03bc" + "M)",
        #yaxis_title="Time (h)",
        #zaxis_title="Normalized cell count",
        xaxis_title="",
        yaxis_title="",
        zaxis_title="",
        xaxis=dict(
            tickmode="array",
            tickvals=doses,
            ticktext=xticks,
            zeroline=False
        ),
        zaxis=dict(
            range=[0, max_z + 0.05],
            tickvals=[""] + list(np.arange((max_z + 0.05) + 1))[1:],
            autorange=False
        )
    ),
    scene_camera={"eye": {"x": 1.8, "y": -1.8, "z": .5}}
)

fig.write_image("../../visualization_figures/vus_norm_visualization.pdf", scale=3)