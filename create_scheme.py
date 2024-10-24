import numpy as np

Vectors = [[-1.000000, 0.000000, 0.000000], 
           [0.000000, 0.707107, 0.707107],
           [0.000000, 0.707107, -0.707107],
           [0.000000, 1.000000, 0.000000 ],
           [-0.707107, 0.000000, -0.707107],
           [0.707107, 0.000000, -0.707107],
           [0.000000, 0.000000, 1.000000],
           [-0.707107, -0.707107, 0.000000],
           [0.707107, -0.707107, 0.000000]] 

bs = [1, 2]
delta = 4 #ms
Deltas = [9.5, 15, 20, 25, 30] #ms
Tex = 54 #ms
TE = 48 #ms

#     x y z G Delta delta TE
# v1 G1
# v1 G2
giro  = 2.6751525e8 * 1e-3 # Gyromagnetic radio [rad/(ms*T)]

f = open("results/T_ex/PGSE_9_dir_3_b.scheme", "w")
f.write("VERSION: STEJSKALTANNER\n")
for vector in Vectors:
    for b in bs:
        for Delta in Deltas:
            G = np.sqrt(b/(pow(giro * delta, 2) * (Delta - delta/3))) # [T/um]
            f.write(f"{vector[0]} {vector[1]} {vector[2]} {G / 1e-6} {delta / 1e3} {Delta / 1e3} {TE}\n")
            b = pow(G * giro * delta, 2) * (Delta - delta/3) # [ms/um²]
            print(b)

f.close
