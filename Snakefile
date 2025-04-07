import os

configfile:
    "configfiles/config.yaml"

data_dir = config["data_dir"]
output_dir = config["output_dir"]
image_file_ext = config["image_file_ext"]
metric = config["metric"]
start_time = config["start_time"]
extrapolate_until_time = config["extrapolate_until_time"]
times = config["times"]
daily = config["daily"]
input_files = [f for f in os.listdir(data_dir) if not f.startswith(".")]

def aggregate_fit_surface_results(wildcards):
    checkpoints.preprocessing.get(**wildcards)
    sample, = glob_wildcards(f"{output_dir}/start_time{start_time}h/separate_files_3D/{{sample}}.csv")
    return expand(f"{output_dir}/start_time{start_time}h/3D/extrapolate_until_{extrapolate_until_time}h/tmp/vus/{{sample}}_daily{{daily}}.csv", sample=sample, daily=daily)

def aggregate_plot_surface(wildcards):
    checkpoints.preprocessing.get(**wildcards)
    sample, = glob_wildcards(f"{output_dir}/start_time{start_time}h/separate_files_3D/{{sample}}.csv")
    return expand(f"{output_dir}/start_time{start_time}h/3D/extrapolate_until_{extrapolate_until_time}h/plots/{{sample}}_time{{time}}_daily{{daily}}.{image_file_ext}", sample=sample, time=times, daily=daily)

rule all:
    input:
        aggregate_plot_surface,
        expand(f"{output_dir}/start_time{start_time}h/3D/extrapolate_until_{extrapolate_until_time}h/fit_surface_results_sorted_by_vus_decrease_after_{{time}}_daily_{{daily}}.csv", time=[t for t in times if t < extrapolate_until_time], daily=daily),
        f"{output_dir}/start_time{start_time}h/3D/extrapolate_until_{extrapolate_until_time}h/fit_surface_mean_scores.txt"

# Preprocessing

checkpoint preprocessing:
    input:
        expand(f"{{data_dir}}/{{input_file}}", data_dir=data_dir, input_file=input_files)
    conda:
        "env.yaml"
    output:
        directory(f"{output_dir}/start_time{start_time}h/separate_files_3D")
    script:
        "scripts/preprocessing.py"

# 3D Surface Rules

rule fit_surface:
    input:
        f"{output_dir}/start_time{start_time}h/separate_files_3D/{{sample}}.csv"
    conda:
        "env.yaml"
    output:
        f"{output_dir}/start_time{start_time}h/3D/extrapolate_until_{extrapolate_until_time}h/tmp/fit_surface/{{sample}}_time{{time}}_daily{{daily}}.csv"
    script:
        "scripts/fit_surface.py"

rule calculate_vus:
    input:
        expand(f"{output_dir}/start_time{start_time}h/3D/extrapolate_until_{extrapolate_until_time}h/tmp/fit_surface/{{sample}}_time{{time}}_daily{{daily}}.csv", time=times, allow_missing=True)
    conda:
        "env.yaml"
    output:
        f"{output_dir}/start_time{start_time}h/3D/extrapolate_until_{extrapolate_until_time}h/tmp/vus/{{sample}}_daily{{daily}}.csv"
    script:
        "scripts/calculate_vus.py"

rule plot_surface:
    input:
        f"{output_dir}/start_time{start_time}h/separate_files_3D/{{sample}}.csv",
        expand(f"{output_dir}/start_time{start_time}h/3D/extrapolate_until_{extrapolate_until_time}h/tmp/fit_surface/{{sample}}_time{{time}}_daily{{daily}}.csv", time=times, daily=daily, allow_missing=True)
    output:
        expand(f"{output_dir}/start_time{start_time}h/3D/extrapolate_until_{extrapolate_until_time}h/plots/{{sample}}_time{{time}}_daily{{daily}}.{image_file_ext}", time=times, daily=daily, allow_missing=True)
    conda:
        "env.yaml"
    script:
        "scripts/plot_surface.py"

rule merge_fit_surface_files:
    input:
        aggregate_fit_surface_results
    conda:
        "env.yaml"
    output:
        f"{output_dir}/start_time{start_time}h/3D/extrapolate_until_{extrapolate_until_time}h/fit_surface_results.csv"
    script:
        "scripts/merge_files.py"

rule sort_by_vus_decrease:
    input:
        f"{output_dir}/start_time{start_time}h/3D/extrapolate_until_{extrapolate_until_time}h/fit_surface_results.csv"
    conda:
        "env.yaml"
    output:
        f"{output_dir}/start_time{start_time}h/3D/extrapolate_until_{extrapolate_until_time}h/fit_surface_results_sorted_by_vus_decrease_after_{{time}}_daily_{{daily}}.csv"
    script:
        "scripts/sort_by_vus_decrease.py"

rule mean_scores_surface:
    input:
        f"{output_dir}/start_time{start_time}h/3D/extrapolate_until_{extrapolate_until_time}h/fit_surface_results.csv"
    conda:
        "env.yaml"
    output:
        f"{output_dir}/start_time{start_time}h/3D/extrapolate_until_{extrapolate_until_time}h/fit_surface_mean_scores.txt"
    script:
        "scripts/mean_scores.py"
