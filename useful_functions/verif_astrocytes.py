import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from DKI import  read_and_extract_dwi, array_to_nifti_replicated, dki_from_file

if __name__ == "__main__":
    path_to_dwi = "/home/localadmin/Documents/MCDS/Permeable_MCDS/output/verif_astrocytes/no_astrocytes_extra_DWI.bfloat"
    scheme_file = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/SMI.scheme"

    _, MD_extra, AD_extra, RD_extra,_,_, _= dki_from_file(path_to_dwi, scheme_file)

    print("MD extra: ", MD_extra)
    print("AD extra: ", AD_extra)
    print("RD extra: ", RD_extra)