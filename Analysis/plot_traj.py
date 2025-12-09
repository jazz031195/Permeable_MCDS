import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

traj = "/home/localadmin/Bureau/_rep_03_0.traj.txt"
t = pd.read_csv(traj, sep=" ")
x = t.x.values
y = t.y.values
z = t.z.values

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(x, y, z)
plt.show()
print(t)