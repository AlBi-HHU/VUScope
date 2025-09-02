import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

import utils

steps_surface = 5
color = px.colors.qualitative.Plotly[0]
colorscale = "Plotly3"
colorscale2 = [[0, "grey"], [1, "grey"]]

separate_file_sample = pd.read_csv(snakemake.input[0], index_col=0)
dfs_to_merge = []
for df in snakemake.input[1:]:
    dfs_to_merge.append(pd.read_csv(df, index_col=0))
df_merged = pd.concat(dfs_to_merge).reset_index(drop=True)
times = df_merged["time"].to_numpy()
max_inhibition_time = max(times) # maximal inhibition time in Incucyte data
start_time = snakemake.config["start_time"]/24
extrapolate_until_time = snakemake.config["extrapolate_until_time"]/24

# Determine maximum height (z-axis value) to use for all plots
max_z = 0
for i in range(df_merged.shape[0]):
    min_dose, max_dose = df_merged.iloc[i, 2:4]
    max_time = extrapolate_until_time - start_time # Assuming 1. extrapolate_until_time >= max_inhibition_time and 2. the highest time always yields a higher cell count than any previous time for at least one of the doses
    best_params = np.array(df_merged.iloc[i, 6:-2])
    dose_steps = np.linspace(min_dose, max_dose, len(separate_file_sample["dose"].unique())*steps_surface + 1)
    time_steps = np.linspace(0, max_time, len(separate_file_sample["time"].unique())*steps_surface + 1)
    norm_cell_count_steps = np.array([utils.dose_time_response_model(best_params, (dose_steps, time_step)) for time_step in time_steps])
    max_z = np.max([max_z, np.max(norm_cell_count_steps)])

# Plot
for i in range(df_merged.shape[0]):
    cell_line, drug = df_merged.iloc[i, :2]
    min_dose, max_dose = df_merged.iloc[i, 2:4]
    time = df_merged.iloc[i, 4]
    daily = df_merged.iloc[i, 5]
    best_params = np.array(df_merged.iloc[i, 6:-2])
    
    fig = go.Figure()
    xticks = 10**separate_file_sample["dose"]
    xticks = [f"{x:.3f}" for x in xticks]
    yticks = separate_file_sample["time"]*24
    yticks = [f"{int(y)}" for y in yticks]
    
    # Draw data points
    max_time = extrapolate_until_time - start_time # Assuming extrapolate_until_time >= max_inhibition_time
    time_corrected = time
    if time == max_inhibition_time:
        time_corrected = time - start_time
    mask = separate_file_sample["time"].isin(list(np.arange(int(time_corrected*24) + 1)/24))
    if daily:
        mask = separate_file_sample["time"].isin(list(np.arange(0, int(time_corrected*24) + 1, 24)/24))
    fig.add_trace(go.Scatter3d(
        x=separate_file_sample["dose"][mask],
        y=separate_file_sample["time"][mask],
        z=separate_file_sample["norm_cell_count"][mask],
        name=f"{cell_line}<br>{drug}",
        marker=dict(size=1, color=color),
        mode="markers",
        showlegend=False
    ))
    
    # Draw surface    
    max_time = extrapolate_until_time - start_time
    dose_steps = np.linspace(min_dose, max_dose, len(separate_file_sample["dose"].unique())*steps_surface + 1)
    time_steps = np.linspace(0, max_time, len(separate_file_sample["time"].unique())*steps_surface + 1)
    norm_cell_count_steps = np.array([utils.dose_time_response_model(best_params, (dose_steps, time_step)) for time_step in time_steps])
    
    if time == max_inhibition_time:
        fig.add_trace(go.Surface(
            x=dose_steps,
            y=time_steps,
            z=norm_cell_count_steps,
            name=f"{cell_line}<br>{drug}",
            opacity=0.5,
            colorscale=colorscale,
            showscale=False
        ))
    else:
        mask1 = time_steps <= time_corrected
        mask2 = time_steps > time_corrected
        fig.add_trace(go.Surface(
            x=dose_steps,
            y=time_steps[mask1],
            z=norm_cell_count_steps[mask1],
            name=f"{cell_line}<br>{drug}",
            opacity=0.5,
            colorscale=colorscale,
            showscale=False
        ))
        
        fig.add_trace(go.Surface(
            x=dose_steps,
            y=time_steps[mask2],
            z=norm_cell_count_steps[mask2],
            name=f"{cell_line}<br>{drug}",
            opacity=0.5,
            colorscale=colorscale2,
            showscale=False
        ))

    fig.update_layout(
        scene=dict(
            xaxis_title="Concentration (" + u"\u03bc" + "M)",
            yaxis_title="Time (h)",
            zaxis_title="Normalized<br> cell count",
            #xaxis_title="",
            #yaxis_title="",
            #zaxis_title="",
            xaxis=dict(
                tickmode="array",
                tickvals=separate_file_sample["dose"],
                ticktext=xticks,
                zeroline=False
            ),
            yaxis=dict(
                tickmode="array",
                tickvals=np.arange(0, int(max_inhibition_time*24) + 1, 24)/24,
                ticktext=np.arange(0, int(max_inhibition_time*24) + 1, 24),
                range=[0, max_inhibition_time]
            ),
            zaxis=dict(
                range=[0, max_z + 0.05],
                tickvals=[""] + list(np.arange((max_z + 0.05) + 1))[1:],
                autorange=False
            ),
            aspectmode="manual",
            aspectratio=dict(x=1, y=1, z=1),
        ),
        scene_camera={"eye": {"x": 1.8, "y": -2.0, "z": 0.4}},
    )
    
    fig.write_image(snakemake.output[i], scale=3)
        