import numpy as np
import matplotlib.pyplot as plt

axis= np.asarray([-1.0, -1.0, 0.0])
origin = np.asarray([0.0, 0.0, 0.0])
plt.plot([axis[0], origin[0]], [axis[1], origin[1]], c="red")

for _ in range(20):

    axis_angle = np.arctan2(axis[1], axis[0])
    temp_angle = np.pi*((1 - 0.95)*np.random.uniform() + (-0.5 + 0.95/2)) + axis_angle;
    new_axis = np.asarray([np.cos(temp_angle), np.sin(temp_angle), 0])

    plt.plot([new_axis[0], origin[0]], [new_axis[1], origin[1]], c="green")

plt.xlim(-2, 2)
plt.ylim(-2, 2)
plt.show()