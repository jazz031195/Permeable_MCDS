import numpy as np
from useful_functions import read_binary_file



if __name__ == "__main__":
    
    file_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/permeable_axons/perm_0_DWI_magnitude_SNR_50.bfloat"

    data = read_binary_file(file_path)

    print(data)
