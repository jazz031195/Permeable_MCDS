import numpy as np

T = 100
dt = 87/T # [ms]
filename = "/Users/ideriedm/Documents/MCDS_perm/Permeable_MCDS/instructions/ISMRM24/Benchmark/overlap4/n2/N_50000_T_10000/_rep_08_0.idx"

data = np.fromfile(filename, dtype="float32")
data = data.reshape((-1, 3))

idx_soma = data[:, 0]
idx_soma[5:8] = 0

changes = np.diff(idx_soma) != 0
change_indices = np.where(changes)[0] + 1

segments = np.split(idx_soma, change_indices)

d = {'soma_dendrites_t': [], 'dendrites_soma_t': []}

for segment in segments:
    print(segment)
    value = segment[0]
    duration = len(segment)
    if value == -1:
        d['dendrites_soma_t'].append(duration*dt)
    else:
        d['soma_dendrites_t'].append(duration*dt)

    print(f"Value {value} lasted {duration} steps")

print(d)

