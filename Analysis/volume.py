import numpy as np

r_soma           = 4e-6 # [m]

nb_dendrites             = 7
length_segment           = 10e-6 # [m]
nb_segments_per_dendrite = 7
radius_dendrite          = 0.8e-6 # [m]

volume_neurites  = np.pi * (radius_dendrite ** 2) * length_segment * nb_segments_per_dendrite * nb_dendrites # in [m³] (3 branching)
volume_soma      = 4/3 * np.pi * r_soma**3 # in [m³]
print("soma volume {:e}".format((volume_soma*1e18)))
print("neurites volume {:e}".format((volume_neurites*1e18)))

volume_soma      = volume_soma * 1e18 # in [um³]
volume_neurites  = volume_neurites * 1e18 # in [um³]
volume_neuron    = volume_neurites + volume_soma
neurite_fraction = volume_neurites / volume_neuron
soma_fraction    = volume_soma / volume_neuron
print("soma volume {:e} in um³".format((volume_soma)))
print("neurites volume {:e}".format((volume_neurites)))
print("neuron {:e}".format((volume_neuron)))

print("\n")
print("soma fraction {:e}".format(soma_fraction))
print("neurites fraction {:e}\n".format(neurite_fraction))

td = 127 # in [ms]
T  = 15500 # Number of timesteps, 15526->0.16um, 6944->0.24um, 3906->0.32um
D0 = 2 # in [um²/ms]
print(f"step length 3D {np.sqrt(2*td/T*D0)} um, step length MCDS {np.sqrt(6*td/T*D0)} um, time step {td/T} s")
print(f"step length 1D {np.sqrt(2*td/T*D0)} um, step length MCDS {np.sqrt(2*td/T*D0)} um, time step {td/T} s")
print(f"Diffusion distance {np.sqrt(2*2*100)} um")
print(f"Diffusion distance {np.sqrt(2*2*60)} um")