import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x = np.linspace(0,1)
y = np.sin(2 * np.pi * x)
for z in np.linspace(0,1,10):
    ax.plot(x, y, z, zdir='x')

ax.view_init(15, 15)

plt.show()
