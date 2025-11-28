import numpy as np
import pandas as pd
from pathlib import Path

input_folder = Path("results/neuron_tortuous_beaded_3/")

traj_files = sorted(input_folder.glob("*traj*.txt"))

dfs = []

for file in traj_files:
    df = pd.read_csv(file, delim_whitespace=True)
    dfs.append(df)

df_all = pd.concat(dfs, ignore_index=True)

df_all = df_all[df_all["step_number"] == 1]
df_traj_init_pos = df_all[["x", "y", "z"]]

df_traj_init_pos = df_traj_init_pos.sample(frac=1).reset_index(drop=True)

df_traj_init_pos.to_csv(input_folder / "init_positions.txt", sep=" ", index=False, header=False)