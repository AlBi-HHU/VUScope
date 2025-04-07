# VUScope
This repository contains the Snakemake pipelines for "VUScope: a mathematical model for evaluating image-based drug response measurements and predicting long-term incubation outcomes".

## Installation
- Clone/Download this repository.
- Install Miniconda.
- Follow the instructions on the [Snakemake installation guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) to install a Snakemake environment in conda. (We used Snakemake version 7.25.0. The versions of the packages are listed in env.yaml.)

## Data
Example data is given in `input`. For the 3D fits, use files formatted like the files in `data` and `protocols` (should be more or less the default output format of Incucyte, but non-Incucyte data could also be formatted like this to be processed). For the 2D fits, use files formatted like the files in `data_celltiterglo`.

## Execution
Activate the Snakemake environment in conda with the command `conda activate snakemake`.

For the main 3D workflow: Prepare the configuration file `configfiles/config.yaml` (see the "Configuration Files" section below). Open the terminal. Navigate to the cloned/downloaded git repo directory. Run:
```bash
snakemake --cores X --use-conda
```
or, if mamba causes problems, add the `--conda-frontend conda` flag:
```bash
snakemake --cores X --use-conda --conda-frontend conda
```
`X` is the number of cores you wish to use. You need to use at least `X=1`. Using multiple cores will parallelize the process. If `X` is greater than the number of cores on your machine, Snakemake will automatically set `X` to the number of cores on your machine. You can set `X=all` to use all your cores. Using the maximum number of cores on your machine means you will likely not be able to use the machine for anything else until the pipeline finished its run. Note: There may appear warnings if overflows occur during calculations, but they can be safely ignored.

For the 2D workflow: Prepare the configuration file `configfiles/config_2D.yaml` (see the "Configuration Files" section below). Open the terminal. Navigate to the cloned/downloaded git repo directory. Run:
```bash
snakemake --cores X --use-conda -s Snakefile_2D
```
and add `--conda-frontend conda` as above if needed.

After running the 3D and 2D workflows, you can get additional result plots and files: Prepare the configuration file `configfiles/config_results.yaml` (see the "Configuration Files" section below). Open the terminal. Navigate to the cloned/downloaded git repo directory. Run:
```bash
snakemake --cores X --use-conda -s Snakefile_results
```
and add `--conda-frontend conda` as above if needed.

## Configuration Files
`configfiles/config.yaml` has the following configurations:
```yaml
data_dir: input/data # directory name with the Incucyte data (or non-Incucyte data brought into this format)
protocol_dir: input/protocols # directory name with the protocol data
output_dir: output # directory name where the output should be stored

image_file_ext: pdf # file extension for plots, can be any file extension supported by plotly
metric: rmse # can be any other metric, but must be implemented in the function "scoring_function(z, z_pred, metric)" in "scripts/utils.py"

start_time: 6 # remove all times that are less than start time, useful if cells have not scattered across the well immediately
extrapolate_until_time: 120 # for any time in the list below named times, predict the results until this value here

times: [72, 120] # fit the data for datapoints up until each time in this list; any time greater than the maximal inhibition time is treated as the maximal inhibition time
daily: [False, True] # whether to fit using all datapoints (False) or only datapoints in 24h-intervals
```

`configfiles/config_2D.yaml` has the following configurations:
```yaml
data_dir: input/data_celltiterglo # directory name with the traditional high-throughput screening data
output_dir: output # directory name where the output should be stored

image_file_ext: pdf # file extension for plots, can be any file extension supported by plotly
metric: rmse # can be any other metric, but must be implemented in the function "scoring_function(z, z_pred, metric)" in "scripts/utils.py"

start_time: 0 # can be ignored, only necessary for using the same scripts of the main 3D workflow
extrapolate_until_time: 120 # for any time in the list below named times, predict the results until this value here

times: [72, 120] # fit the data for datapoints at each time in this list; any time greater than the maximal inhibition time is treated as the maximal inhibition time
```

`configfiles/config_results.yaml` has the following configurations:
```yaml
output_dir: output # directory name where the output of the 3D and 2D workflows are stored

image_file_ext: pdf # file extension for plots used by the 3D and 2D workflows

start_time_incucyte: 6 # start time used in "configfiles/config.yaml"
start_time_celltiterglo: 0 # start time used in "configfiles/config_2D.yaml"

intermediate_inhibition_time: 72 # a time in times that appears in both "configfiles/config.yaml" and "configfiles/config_2D.yaml" that is less than the maximal inhibition time
max_inhibition_time: 120 # the maximal inhibition time

daily: [False, True] # whether to get result files for output that was fit using all datapoints (False) or only datapoints in 24h-intervals

grayscale: True # True: gray plots without background, False: colored plots with background
```

## Documentation
Here is a short description of the scripts in `scripts`.

Main 3D workflow:
- `preprocessing.py`: Merges the Incucyte output of a cell line and its respective protocol. Removes all values less than the start time given in the configuration file. Normalizes the cell counts and log10-transfroms them.
- `fit_surface.py`: Fits the dose-time-response surface.
- `calculate_vus.py`: Calculates the VUS by averaging the numerical VUS approximation above and below the surface.
- `plot_surface.py`: Plots the dose-time-response surface.
- `merge_files.py`: Merges files that were separate for each cell line, drug, and time into one.
- `sort_by_vus_decrease.py`: Sorts the results table such that the most promising drugs that show a stronger effect with more incubation time are at the top of the table.
- `mean_scores.py`: Shows the metric of the fit specified in the configuration file.

2D workflow:
- `scripts_2D/preprocessing_celltiterglo.py`: Brings the CellTiter-Glo output into the same format as the Incucyte data after `preprocessing.py`.
- `scripts_2D/fit_curve.py`: Fits the dose-response curve.
- `scripts_2D/calculate_auc.py`: Calculates the AUC analytically.
- `scripts_2D/plot_curve.py`: Plots the dose-response curve.
- `merge_files.py`: See main 3D workflow.

Results workflow:
- `scripts_results/norm_cell_count_distribution.py`: Plots the distribution of normalized cell counts. This plot is not in the manuscript and only helps to see the range of values to get a feeling whether a certain metric value is good or bad (e.g., whether an RMSE of 0.4 is actually good (if the range of values is large, for example between 0 and 15) or bad (if the range of values is between 0 and 1)).
- `scripts_results/trend_plot.py`: Plots the trend whether more incubation time is predicted to improve drug efficacy or actually improves drug efficacy, see the manuscript.
- `scripts_results/correlations_and_mapes_vus.py`: Plots and lists the Pearson and concordance correlation coefficients and mean absolute percentage errors of the VUS values.
- `scripts_results/correlations_and_mapes_auc.py`: Plots and lists the Pearson and concordance correlation coefficients and mean absolute percentage errors of the AUC values.
- `scripts_results/merge_txt_files.py`: Similar to `merge_files.py`.
