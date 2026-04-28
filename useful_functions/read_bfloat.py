import numpy as np
from useful_functions import read_binary_file



if __name__ == "__main__":
    
    file_path = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/SMI_pred/straight_cyl/combined/f_0.2/DWI_total.bfloat"

    data = read_binary_file(file_path)

    print(data)
