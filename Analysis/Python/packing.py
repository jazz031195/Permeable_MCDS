"""
    File to compare the mesh with the analytical solutions for different 
    decimation. Plots the mean diffusivity (MD) for different decimations
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
if sys.platform == "linux":
    sys.path.insert(1, '/home/localadmin/Documents/analytical_formula/')
else:
    sys.path.insert(1, '/Users/ideriedm/Documents/analytical_formula/')

import pandas as pd
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")
from utils import get_bvals, get_bvectors, calculate_DKI, get_dwi, get_psge, create_data, analytical_solutions


scheme_file_path = "docs/scheme_files/PGSE_21_dir_12_b.scheme"
filename     = '_rep_12_DWI_img.txt'
# Name of the experience
name         = ('_').join(filename.split('_')[:-1])
extension    = filename.split('_')[-1].split('.')[-1]
SNR          = np.inf
data_one_exp = create_data(Path('/Users/ideriedm/Documents/MCDS_perm/Permeable_MCDS/results/packing/'), 
                           SNR, 
                           name, 
                           extension, 
                           scheme_file_path)

print(data_one_exp.columns)
print(data_one_exp[['adc [ms/um²]', 'FA', 'MD', 'AD', 'RD', 'b [ms/um²]']])