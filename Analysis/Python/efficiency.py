import matplotlib.pyplot as plt
import numpy as np

sizes = [6.6e-3, 4.6, 9.3, 24.3, 49.2, 74.2, 387.3]
plt.bar(list(range(len(sizes))), sizes)
plt.ylabel("Neuron file size [Mo]")
plt.xticks(list(range(len(sizes))), ["Spheres", "Mesh 0.05", "Mesh 0.10", "Mesh 0.25", "Mesh 0.50", "Mesh 0.75", "Mesh 1"])
plt.show()