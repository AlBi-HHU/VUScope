import os

import numpy as np
import pandas as pd

def clean_well_names(well_list):
    """
    Standardizes well names by converting formats like A01 to A1.
    """
    for i in range(len(well_list)):
        letter, count = well_list[i][0], str(int(well_list[i][1:]))
        well_list[i] = letter + count
    return well_list

output_dir = snakemake.config["output_dir"]
protocol_dir = snakemake.config["protocol_dir"]
start_time = int(snakemake.config["start_time"])

dfs_to_merge = []
for input_file in snakemake.input:
    input_parts = input_file.split("/")[-1].split("_")
    cell_line = input_parts[0]
    protocol_name = input_parts[1].split(".")[0]

    # Read and clean data
    header_index = 0
    with open(input_file) as file:
        for line in file:
            if "Elapsed" in line:
                break
            if line.strip():
                header_index += 1

    df_data = pd.read_csv(input_file, delimiter="\t", header=header_index, index_col="Elapsed")
    df_data = df_data.drop("Date Time", axis=1).T
    
    # Remove all times < start_time
    old_time_cols = df_data.columns.to_numpy()
    num_cols_to_drop = np.sum(old_time_cols < start_time)
    new_time_cols = {old_time_cols[i]: old_time_cols[i - num_cols_to_drop] for i in range(num_cols_to_drop, len(old_time_cols))}
    while num_cols_to_drop > 0:
        df_data = df_data.drop(df_data.columns[0], axis=1)
        num_cols_to_drop -= 1
    df_data = df_data.rename(columns=new_time_cols)

    # Normalize cell counts by division by cell count at time 0
    for col in df_data.columns.to_numpy()[1:]:
        df_data[col] /= df_data[0]
    df_data[0] /= df_data[0]

    df_data["CellLine"] = cell_line
    
    # Read protocols
    protocol_file = f"{protocol_dir}/{protocol_name}_tabular_detail_anonymized.xlsx" # for non-anonymized IncuCyte excel sheets, remove "_tabular_detail_anonymized"
    df_protocol = pd.read_excel(protocol_file)
    df_protocol.columns = [col.replace("\n", " ") for col in df_protocol.columns]

    # Clean protocols
    df_protocol["Dispensed well"] = clean_well_names(df_protocol["Dispensed well"].values)
    df_protocol = df_protocol[(df_protocol["Fluid name"] != "DMSO normalization") | (~df_protocol["Dispensed well"].duplicated())]
    df_protocol = df_protocol[["Fluid name", "Dispensed concentration", "Dispensed well"]]
    df_protocol["Protocol"] = protocol_name

    # Merge data and protocol
    merge_columns = ["Fluid name", "CellLine", "Protocol", "Dispensed concentration", "Dispensed well"]
    df = df_data.merge(df_protocol, right_on="Dispensed well", left_index=True)
    df = df.melt(id_vars=merge_columns, var_name="time", value_name="norm_cell_count")

    # log10-transform doses
    df["Dispensed concentration"] = np.log10(df["Dispensed concentration"])
    
    dfs_to_merge.append(df)

# Save all samples (cell-line--drug--protocol triplets)
df_merged = pd.concat(dfs_to_merge)
df_merged = df_merged.rename(columns={
    "Fluid name": "drug",
    "CellLine": "cell_line",
    "Dispensed concentration": "dose"
})
df_merged = df_merged[["cell_line", "drug", "dose", "time", "norm_cell_count"]]
df_merged = df_merged.sort_values(by=["cell_line", "drug", "dose", "time"])

output_path = f"{output_dir}/start_time{start_time}h/separate_files"
if not os.path.exists(output_path):
    os.makedirs(output_path)

all_cell_lines = df_merged["cell_line"].unique()
for cell_line in all_cell_lines:
    sample_data = df_merged[df_merged["cell_line"] == cell_line]
    all_drugs = [drug for drug in sample_data["drug"].unique() if drug != "DMSO normalization"]
    for drug in all_drugs:
        sample = sample_data[sample_data["drug"] == drug]
        sample = sample.reset_index(drop=True)
        sample.to_csv(f"{output_path}/{cell_line}_{drug}.csv")