import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats
import random

def shrink_swc(path, new_path, shrink_perc):
    # Open the input file in read mode and the output file in write mode
    with open(path, 'r') as infile, open(new_path, 'w') as outfile:
        lines = infile.readlines()
        for i in range(2):
            outfile.write(lines[i])
        lines = lines[2:]
        for i in range(len(lines)):
            coords = lines[i].split(' ')
            if len(coords) > 2 and len(coords) < 6:
                x = coords[0] 
                y = coords[1] 
                z = coords[2]
                r = float(coords[3])
                new_r = shrink_radius(shrink_perc, r)
                outfile.write(f'{x} {y} {z} {new_r}\n')
            else:
                outfile.write(lines[i])


    
def shrink_radius(perc, swollen_radius):
    return swollen_radius/np.sqrt(1+perc)




if __name__ == "__main__":

    # shrink by 1%
    path = "results/Simu/neurons_list_49_perc.txt"
    new_path = "results/Simu/neurons_list_49_perc_1_perc_shrink.txt"
    shrink_perc = 0.01
    new_lines = shrink_swc(path, new_path, shrink_perc)

