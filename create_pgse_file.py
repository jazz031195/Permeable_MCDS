import itertools
import math

giro = 267.51525e3 #rad/msXT

def generate_combinations_b_vals(x_vals, y_vals, z_vals, b_vals, delta_vals, te_vals, delta_small_vals, output_file):
    with open(output_file, 'w') as f:
        f.write("VERSION: STEJSKALTANNER\n")
        for (x, y, z), b, (Delta, TE, delta) in itertools.product(
                zip(x_vals, y_vals, z_vals), b_vals, zip(delta_vals, te_vals, delta_small_vals)):
            G = (math.sqrt(b/(Delta- (delta/3))))/(delta*giro) #T/m
            # round G 
            G = round(G, 8)
            print("b : ", b)
            print("G : ",G)
            Delta = round(Delta, 8)
            TE = round(TE, 8)
            delta = round(delta, 8)

            f.write(f"{x} {y} {z} {G} {Delta} {delta} {TE}\n")

def generate_combinations_G(x_vals, y_vals, z_vals, g_vals, delta_vals, te_vals, delta_small_vals, output_file):
    with open(output_file, 'w') as f:
        f.write("VERSION: STEJSKALTANNER\n")
        for (x, y, z), G, (Delta, TE, delta) in itertools.product(
                zip(x_vals, y_vals, z_vals), g_vals, zip(delta_vals, te_vals, delta_small_vals)):

            # round G 
            G = round(G, 4)
            Delta = round(Delta, 4)
            TE = round(TE, 4)
            delta = round(delta, 4)

            f.write(f"{x} {y} {z} {G} {Delta} {delta} {TE}\n")

x_vals = [0]
y_vals = [0]
z_vals = [1]
bs = [0, 200, 500, 1000, 2000, 5000] # s/mm²
print(bs)
Deltas = [0.1013, 0.0463, 0.0263, 0.0173, 0.0123]
deltas = [0.004]*len(Deltas)
tes = [Deltas[i]+deltas[i]+0.005 for i in range(len(Deltas))]

output_file = "/home/localadmin/Documents/MCDS/Permeable_MCDS/instructions/scheme/time_dependence.scheme"
generate_combinations_b_vals(x_vals, y_vals, z_vals, bs, Deltas, tes, deltas, output_file)

