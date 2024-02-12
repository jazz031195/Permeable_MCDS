# To compare ADC between multiple DWI files with varying N values

import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

cur_path = os.getcwd()
giro = 2.6751525e8 #Gyromagnetic radio given in rad/(ms*T)
scheme_file = cur_path + "/MCDC_Simulator_public-master/docs/scheme_files/PGSE_sample_scheme_new.scheme"

def get_dwi_array(dwi_path):
    # create an array with dwi values
    signal = []
    with open(dwi_path) as f:
        # read each line in the file
        # convert the line to a float and append it to the signal array
        [signal.append(float(line)) for line in f.readlines()]
    return np.array(signal)

def get_psge_data():
    # create an empty DataFrame with specified column names
    data_dwi = pd.DataFrame(columns=["x", "y", "z", "G", "Delta", "delta", "TE"])
    
    # create empty lists for each column
    x, y, z, G, Delta, delta, TE = [], [], [], [], [], [], []
    
    with open(scheme_file) as f:
        for line in f.readlines():
            if len(line.split(' ')) > 2:
                for e, element in enumerate(line.split(' ')):
                    if e == 0:
                        x.append(float(element))  # add x-coordinate value to the x list
                    elif e == 1:
                        y.append(float(element))  # add y-coordinate value to the y list
                    elif e == 2:
                        z.append(float(element))  # add z-coordinate value to the z list
                    elif e == 3:
                        G.append(float(element) * 1e-3)  # add G value (multiplied by 1e-3) to the G list
                    elif e == 4:
                        Delta.append(float(element))  # add Delta value to the Delta list
                    elif e == 5:
                        delta.append(float(element))  # add delta value to the delta list
                    elif e == 6:
                        TE.append(float(element[:-1]))  # add TE value (excluding the last character) to the TE list
    
    # assign the lists to respective columns in the DataFrame
    data_dwi["x"] = x
    data_dwi["y"] = y
    data_dwi["z"] = z
    data_dwi["G"] = G
    data_dwi["Delta"] = Delta
    data_dwi["delta"] = delta
    data_dwi["TE"] = TE
    
    # calculate a new column "b [ms/um²]" based on existing columns
    data_dwi["b [ms/um²]"] = pow(data_dwi["G"] * giro * data_dwi["delta"], 2) * \
                             (data_dwi["Delta"] - data_dwi["delta"] / 3) / 1000
    
    return data_dwi

def create_data(dwi_path, b0):
    # Get the DWI signal array from the given file path
    dwi_signal = get_dwi_array(dwi_path)
    
    # Get the psge data DataFrame
    data_psge = get_psge_data()
    
    # Add the DWI signal array as a new column named "DWI" in the data_psge DataFrame
    data_psge["DWI"] = list(dwi_signal)

    # Filter data based on x, y, and z coordinates greater than 0.0
    data_x = data_psge.loc[data_psge['x'] > 0.0]
    data_y = data_psge.loc[data_psge['y'] > 0.0]
    data_z = data_psge.loc[data_psge['z'] > 0.0]

    # Create a list of filtered dataframes
    datas = [data_x, data_y, data_z]

    for i in range(len(datas)):
        # Round b values to 1 decimal place
        datas[i]["b [ms/um²]"] = [round(n, 1) for n in list(datas[i]["b [ms/um²]"])]

        # Get the DWI(b0) value for the given b0
        Sb0 = list(datas[i].loc[datas[i]["b [ms/um²]"] == b0]["DWI"])[0]

        # Calculate log(Sb/So) for each DWI value in the dataframe
        signal = list(map(lambda Sb: np.log(Sb / Sb0), list(datas[i]["DWI"])))

        # Calculate adc [um²/ms] for each DWI value in the dataframe
        adc = list(map(lambda b, Sb: -np.log(Sb / Sb0) / (b - b0) if b != b0 else np.nan,
                       list(datas[i]["b [ms/um²]"]), list(datas[i]["DWI"])))

        # Add the calculated columns to the dataframe
        datas[i]["log(Sb/So)"] = signal
        datas[i]["adc [um²/ms]"] = adc

    # Concatenate the filtered dataframes into a single dataframe
    data_dwi = pd.concat(datas)

    return data_dwi

def get_adc(dwi_path, b, b0):
    data = create_data(dwi_path, b0)
    axes = ["x","y","z"]
    orientations = ["radial","radial","axial"]
    adcs = []
    for a in axes:
        ax_data = data.loc[data[a]>0]
        ax_data = ax_data.loc[ax_data["b [ms/um²]"] == b]
        adcs.append(list(ax_data["adc [um²/ms]"])[0])
    new_data = pd.DataFrame(columns = ["axis", "adc [um²/ms]"])
    new_data["orientations"] = orientations
    new_data["axis"] = axes
    new_data["adc [um²/ms]"] = adcs
    return new_data


def get_axons_array(file):
    axons = []  # List to store axons
    spheres = []  # List to store spheres

    with open(file) as f:
        for line in f.readlines():
            if len(line.split(' ')) > 4:
                # If the line has more than 4 elements, it represents a sphere
                # Convert the elements to float and add them to the spheres list
                spheres.append([float(i) for i in line.split(' ')[:]])
            else:
                # If the line does not have more than 4 elements, it represents the end of a set of spheres
                # Check if spheres list is not empty (to avoid adding empty lists)
                if len(spheres) != 0:
                    axons.append(spheres)  # Add the spheres list to the axons list
                spheres = []  # Reset the spheres list for the next set

    return axons

def plot_radii(file):
    # Get the axons array from the given file
    axons = get_axons_array(file)

    for i in range(len(axons)):
        radius = []  # List to store radii
        spheres = axons[i]  # Get the spheres for the current axon
        for ii in range(len(spheres)):
            radius.append(spheres[ii][3] * 1000)  # Extract and convert the radius (index 3) to micrometers
        if len(list(dict.fromkeys(radius))) > 200:
            # If the number of unique radii is greater than 200, plot the radii
            plt.plot(radius)

    plt.ylabel("[um]")  # Set the y-axis label as "[um]"
    plt.show()  # Display the plot

def plot_tortuosity_(file):
    axons = get_axons_array(file)
    tortuosites =[] 
    radii =[] 
    for axon in axons :
        length = 0 
        mean_radius = 0
        for i, sphere in enumerate(axon):
            if i>0:
                length += np.linalg.norm(np.array([axon[i-1][0], axon[i-1][1], axon[i-1][2]])- np.array([axon[i][0], axon[i][1], axon[i][2]]))
            mean_radius += axon[i][3]*1000
        dist = np.linalg.norm(np.array([axon[0][0],axon[0][1] ,axon[0][2]])- np.array([axon[-1][0],axon[-1][1] ,axon[-1][2]]))
        tortuosites.append(length/dist)
        radii.append(mean_radius/len(axon))
    df = pd.DataFrame()
    df["radius [um]"]  = radii 
    df["tortuosity"]  = tortuosites
    g = sns.JointGrid(data=df, x="radius [um]", y="tortuosity", space=0)
    g.plot_joint(sns.kdeplot,
             fill=True,
              cmap="rocket")
    g.plot_marginals(sns.histplot, color="#03051A", alpha=1, bins=25)
    plt.show()


def plot_tortuosity(file):
    axons = get_axons_array(file)
    dfs =[] 
    for j in range(int(len(axons)/100)):
        i = j*50
        xs = [] 
        ys =[]  
        spheres = axons[i]
        for ii in range(1,len(spheres)):
            xs.append(spheres[ii][0]-spheres[ii-1][0])
            ys.append(spheres[ii][1]- spheres[ii-1][1])
        df = pd.DataFrame()
        df["x[um]"]=xs
        df["y[um]"]=ys
        df["index"] = [i]*len(xs)   
        dfs.append(df)
    dfs_ = pd.concat(dfs)
    dfs_ = dfs_.reset_index()
    g = sns.JointGrid(data=dfs_, x="x[um]", y="y[um]")
    g.plot_joint(sns.kdeplot,
             fill=True, 
             thresh=0, levels=100, cmap="rocket")
    g.plot_marginals(sns.histplot, color="#03051A", alpha=1, bins=25)
    plt.show()

def plot_tortuosity_angle(file):
    axons = get_axons_array(file)
    dfs =[] 
    for j in range(int(len(axons)/100)):
        i = j*50
        xs = [] 
        ys =[]  
        zs =[] 
        spheres = axons[i]
        for ii in range(1,len(spheres)):
            xs.append(spheres[ii][0]-spheres[ii-1][0])
            ys.append(spheres[ii][1]- spheres[ii-1][1])
            zs.append(spheres[ii][2]- spheres[ii-1][2])

        df = pd.DataFrame()
        df["index"] = [i]*len(xs)
        df["theta"]= list(map(lambda x,y:np.arctan2(y,x)/np.pi if np.arctan2(y,x)/np.pi>0 else np.arctan2(y,x)/np.pi+2, xs, ys))
        df["phi"]= list(map(lambda x,y,z:np.arctan(np.sqrt(x**2+y**2)/z)/np.pi if z!= 0 else np.pi/2, xs, ys, zs))  
        dfs.append(df)
        del df
    dfs_ = pd.concat(dfs)
    dfs_ = dfs_.reset_index()
    g = sns.JointGrid(data=dfs_, x="theta", y="phi")
    g.plot_joint(sns.kdeplot,
             fill=True, 
             thresh=0, levels=100, cmap="rocket")
    g.plot_marginals(sns.histplot, color="#03051A", alpha=1, bins=25)
    plt.show()
    

def plot_DWI(intra_rest_path,extra_rest_path, extra_active_path, intra_active_path):
    data_intra_rest = create_data(intra_rest_path)
    data_extra_active = create_data(extra_active_path)
    data_intra_active = create_data(intra_active_path)
    data_extra_rest = create_data(extra_rest_path)
    data_rest = data_intra_rest.add(data_extra_rest)
    data_active = data_intra_active.add(data_extra_active)
    ax = sns.lineplot(data=data_rest, x = "b [ms/um²]", y= "DWI")

    ax = sns.lineplot(data=data_active, x = "b [ms/um²]", y= "DWI")
    plt.show()


def create_rel_df(data):
    data = data.groupby(by = ["swelling percentage [%]","location"]).mean().reset_index()
    data = data.groupby(by = ["swelling percentage [%]"]).sum().reset_index()
    data_0 = data.groupby(by = ["swelling percentage [%]"]).mean().reset_index()
    data_0 = data_0.loc[data_0["swelling percentage [%]"]==0]
    data["Relative_adc"] = list(map(lambda a : a/list(data_0["adc_weighted [um²/ms]"])[0] , list(data["adc_weighted [um²/ms]"])))
    return data

def create_rel_df_(data):
    data = data.groupby(by = ["swelling percentage [%]"]).mean().reset_index()
    data_0 = data.loc[data["swelling percentage [%]"]==0]
    data["Relative_adc"] = list(map(lambda a : a/list(data_0["adc_weighted [um²/ms]"])[0] , list(data["adc_weighted [um²/ms]"])))
    return data

def create_df(straight,no_axial,no_radial, type, icvf_, b1, b0):
        directory = "/home/localadmin/Documents/permeable_MCDS/MCDC_Simulator_public/instructions/axons_vs_cylinders/data"
        locs = ["extra"]
        swell = [0, 30, 50, 70, 100]
        vox_size = 50
        datas = []
        for n in range(2):
            for s in swell:
                for l in locs:
                    if not straight:
                        if type == "axons":
                            if n == 0:
                                file = f"{directory}/axons_vs_cylinders/axons_icvf_{icvf_}_cap_24_vox_50_straight_factor_2_extra_DWI.bfloat"
                            elif n<=10:
                                file = f"{directory}/axons_vs_cylinders/icvf_{icvf_}_vox_{vox_size}_swell_{s}_{l}_rep_0{n-1}_DWI.txt"
                        
                        elif type == "cylinders":
                            if n == 0:
                                file = f"/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/cylinders/icvf_{icvf_}/icvf_{icvf_}_vox_{vox_size}_swell_{s}_{l}_DWI.txt"
                            elif n<=10:
                                file = f"/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/cylinders/icvf_{icvf_}/icvf_{icvf_}_vox_{vox_size}_swell_{s}_{l}_rep_0{n-1}_DWI.txt"
                            else:
                                file = f"/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/cylinders/icvf_{icvf_}/icvf_{icvf_}_vox_{vox_size}_swell_{s}_{l}_rep_{n-1}_DWI.txt"
                        else:
                            print(f"Wrong type : {type}")
                    """
                    else:
                        if type == "axons":
                            if n == 0:
                                file = cur_path + f"/MCDC_Simulator_public-master/instructions/demos/output/{type}/Runs/Straight/icvf_{icvf_}_swell_{s}_straight_{l}_DWI.txt"
                            elif n<=10:
                                file = cur_path + f"/MCDC_Simulator_public-master/instructions/demos/output/{type}/Runs/Straight/icvf_{icvf_}_swell_{s}_straight_{l}_rep_0{n-1}_DWI.txt"
                            else:
                                file = cur_path + f"/MCDC_Simulator_public-master/instructions/demos/output/{type}/Runs/Straight/icvf_{icvf_}_swell_{s}_straight_{l}_rep_{n-1}_DWI.txt"
                        elif type == "cylinders":
                            if n == 0:
                                file = cur_path + f"/MCDC_Simulator_public-master/instructions/demos/output/{type}/Runs/icvf_{icvf_}_swell_{s}_straight_{l}_DWI.txt"
                            elif n<=10:
                                file = cur_path + f"/MCDC_Simulator_public-master/instructions/demos/output/{type}/Runs/icvf_{icvf_}_swell_{s}_straight_{l}_rep_0{n-1}_DWI.txt"
                            else:
                                file = cur_path + f"/MCDC_Simulator_public-master/instructions/demos/output/{type}/Runs/icvf_{icvf_}_swell_{s}_straight_{l}_rep_{n-1}_DWI.txt"
                        else:
                            print(f"Wrong type : {type}")
                    """

                    if (Path(file)).exists():

                        data = get_adc(file, b1, b0)
                        data["repetition"] = [n]*len(data)
                        data["swelling percentage [%]"] = [s]*len(data)
                        data["location"] = [l]*len(data)
                        data["icvf"] = icvf_
    
                        datas.append(data)

                    else:
                        print(f"Careful! File {file} does not exist !")


        data = pd.concat(datas)

        if no_axial :
            data = data.loc[data["orientations"] != "axial"]
        elif no_radial:
             data = data.loc[data["orientations"] != "radial"]

        data = data.reset_index()
 
        # add weighted adc (taking into account ICVF)
        #data["adc_weighted [um²/ms]"] = list(map(lambda x,y,i : x*i if y == "intra" else x*(1-i) , list(data["adc [um²/ms]"] ), list(data["location"]), list(data["icvf"])))
    
        return data

def plot_increase(no_axial, no_radial):
    
    b0 = 0.2
    b1 = 1
    icvf_='0.50'
    """
    type = "cylinders"
    straight = False


    data_cyl = create_df(straight,no_axial,no_radial, type, icvf_, b1, b0)

    type = "axons"
    straight = True
    data_ax_str = create_df(straight,no_axial,no_radial, type, icvf_, b1, b0)
    """
    type = "axons"
    straight = False
    data_ax = create_df(straight,no_axial,no_radial, type, icvf_, b1, b0)

    fig, ax = plt.subplots(figsize=(8, 6))

    data = data_ax
    print(data)
    #data = create_rel_df(data)

    sns.regplot(data = data.loc[data["location"]=="intra"], x = "swelling percentage [%]",y = "adc [um²/ms]",  label = "Tortuous overlapping spheres intra")
    plt.show()
    sns.regplot(data = data.loc[data["location"]=="extra"], x = "swelling percentage [%]",y = "adc [um²/ms]", label = "Tortuous overlapping spheres extra")
    plt.show()
    """
    data = data_cyl
    data = create_rel_df(data)
    sns.regplot(data = data, x = "swelling percentage [%]",y = "Relative_adc", ax= ax, label = "Straight cylinder", ci = None)
    
    data = data_ax_str
    data = create_rel_df(data)
    sns.regplot(data = data, x = "swelling percentage [%]",y = "Relative_adc", ax= ax, label = "Straight overlapping spheres", ci = None)
    # Create legend
    ax.legend()
    if no_axial:
        plt.title("Radial ADC")
    else:
        plt.title("Mean ADC")
    plt.show()


    # Create the subplots
    palette =sns.color_palette("Blues", n_colors=3)
    sns.set_palette(palette)
    fig, ax = plt.subplots(1, 1, figsize=(10, 5))

    data = data_ax
    data_intra = data.loc[data["location"]=="intra"]
    data_extra = data.loc[data["location"]=="extra"]
    data_intra = create_rel_df_(data_intra)
    data_extra = create_rel_df_(data_extra)
    data_tot = create_rel_df(data)


    # Plot the boxplots on the subplots
    sns.lineplot(data = data_intra, x = "swelling percentage [%]",y = "Relative_adc", label = 'Intracellular water molecules', ax= ax)
    sns.lineplot(data = data_extra, x = "swelling percentage [%]",y = "Relative_adc",  label = 'Extracellular water molecules', ax= ax)
    sns.lineplot(data = data_tot, x = "swelling percentage [%]",y = "Relative_adc",  label = 'All water molecules', ax= ax)
    
    ax.legend()

    if no_axial:
        plt.title("RD in Tortuous Axons")
    else:
        plt.title("MD in Tortuous Axons")
    plt.show()
    """

    

    

def get_icvf_cylinders(file):
    area_intra = 0
    with open(file) as f:
        for e,line in enumerate(f.readlines()):
            if e == 5:
                voxel_size = float(line)
            if len(line.split(" "))> 5:
                r = float(line.split(" ")[-2])
                s = float(line.split(" ")[-1])
                if s >0 :
                    R = r*np.sqrt(1.01)
                    area_intra += np.pi*R*R
                else:  
                    area_intra += np.pi*r*r
    total_area = voxel_size*voxel_size
    icvf = area_intra/total_area
    return icvf

def get_icvf_axons(file):
    area_intra = 0
    new_axon = False
    with open(file) as f:
        for e,line in enumerate(f.readlines()):
            if e == 5:
                voxel_size = float(line)
            if len(line.split(" "))> 3 and new_axon:
                r = float(line.split(" ")[-2])
                s = float(line.split(" ")[-1])
                if s >0 :
                    R = r*np.sqrt(1.01)
                    area_intra += np.pi*R*R
                else:  
                    area_intra += np.pi*r*r
                new_axon = False
            if len(line.split(" ")) == 3:
                new_axon = True
    total_area = voxel_size*voxel_size
    icvf = area_intra/total_area
    return icvf

plot_increase(False, False)

        