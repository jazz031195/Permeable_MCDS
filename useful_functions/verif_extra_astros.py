import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from DKI import  dki_from_file


if __name__ == "__main__":
    scheme_file = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/SMI.scheme"
    file_extra = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/verif_astrocytes/astrocytes_0.06_extra_DWI.bfloat"
    # check if the file exists
    if not os.path.isfile(file_extra):
        print(f"File not found: {file_extra}")
    else:
        print(f"File found: {file_extra}")
    _, MD_extra, AD_extra, RD_extra,_,_, _= dki_from_file(file_extra, scheme_file, binary= True)
    print("MD extra: ", MD_extra)
    print("AD extra: ", AD_extra)
    print("RD extra: ", RD_extra)