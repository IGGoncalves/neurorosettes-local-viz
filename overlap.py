import numpy as np
import matplotlib.pyplot as plt
from neurorosettes.physics import get_cylinder_intersection

n1_proximal = np.asarray([-1.0, 2.0, 0.0])
n1_distal = np.asarray([0.0, -2.0, 0.0])
n2_proximal = np.asarray([-3.0, 0.0, 0.0])
n2_distal = np.asarray([3.0, 1.0, 0.0])

def ccw(A,B,C):
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

# Return true if line segments AB and CD intersect
def intersect(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

neurite_axis = np.subtract(n1_distal, n1_proximal)

print(intersect(n1_proximal, n1_distal, n2_proximal, n2_distal))

good_point = get_cylinder_intersection(n1_proximal, n1_distal, n2_proximal, n2_distal)[0]
print(good_point)

fraction = np.linalg.norm(np.subtract(good_point, n1_proximal)) / np.linalg.norm(neurite_axis)

fraction -= 0.05
good_point = n1_proximal + fraction*neurite_axis

plt.plot([n1_proximal[0], good_point[0]], [n1_proximal[1], good_point[1]], c="blue")
plt.plot([n2_proximal[0], n2_distal[0]], [n2_proximal[1], n2_distal[1]], c="red")
plt.show()
