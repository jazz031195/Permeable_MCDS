import numpy as np

for i in range(1, 6):
    folder = f"/home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/branches_single_neuron/straight/overlap4/dendrites/n{i}/" 
    for img in ["_DWI.bfloat", "_rep_00_DWI.bfloat", "_rep_01_DWI.bfloat", "_rep_02_DWI.bfloat", "_rep_03_DWI.bfloat"]:
        print(folder + img)
        print(np.fromfile(folder + img, dtype="float32").shape)