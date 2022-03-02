import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.transforms as transforms

from neurorosettes import physics
from neurorosettes.subcellular import ObjectFactory

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

names = ["Oblique cylinders", "Perpendicular cylinders", "Parallel cylinders"]
springs = [{"coordinates": np.array([-1.0, -3.0, 0.0]), "growth": np.array([1.0, 1.0, 0]), "angle": -45},
           {"coordinates": np.array([-3.0, 0.0, 0.0]), "growth": np.array([1.0, 0.0, 0]), "angle": -90},
           {"coordinates": np.array([2.0, -3.0, 0.0]), "growth": np.array([0.0, 1.0, 0]), "angle": 0}]

for i, ax in enumerate(axes):
    spring1 = ObjectFactory().get_neurite(proximal_position=np.array([-5.0, -5.0, 0.0]),
                                          axis=np.array([0.0, 1.0, 0.0]))

    spring2 = ObjectFactory().get_neurite(proximal_position=springs[i]["coordinates"],
                                          axis=springs[i]["growth"])

    point1, point2 = physics.get_cylinder_intersection(spring1.proximal_point,
                                                       spring1.distal_point,
                                                       spring2.proximal_point,
                                                       spring2.distal_point)

    ax.add_patch(patches.Rectangle((spring1.proximal_point[0] - spring1.mechanics.radius, spring1.proximal_point[1]),
                                   width=spring1.mechanics.radius * 2,
                                   height=spring1.current_length,
                                   edgecolor="firebrick",
                                   facecolor="salmon"))

    rect2 = patches.Rectangle((spring2.proximal_point[0] - spring2.mechanics.radius, spring2.proximal_point[1]),
                              width=spring2.mechanics.radius * 2,
                              height=spring2.current_length,
                              edgecolor="darkkhaki",
                              facecolor="moccasin")

    t2 = transforms.Affine2D().rotate_deg_around(spring2.proximal_point[0],
                                                 spring2.proximal_point[1],
                                                 springs[i]["angle"]) + ax.transData
    rect2.set_transform(t2)
    ax.add_patch(rect2)

    ax.plot([spring1.proximal_point[0], spring1.distal_point[0]],
            [spring1.proximal_point[1], spring1.distal_point[1]],
            color="firebrick", linestyle="--", zorder=1)

    ax.plot([spring2.proximal_point[0], spring2.distal_point[0]],
            [spring2.proximal_point[1], spring2.distal_point[1]],
            color="darkkhaki", linestyle="--", zorder=1)

    ax.add_patch(patches.Circle((point1[0], point1[1]), radius=0.2,
                                fill=True, alpha=1, facecolor="black"))

    ax.add_patch(patches.Circle((point2[0], point2[1]), radius=0.2,
                                fill=True, alpha=1, facecolor="black"))

    ax.scatter([], [], s=10, marker='o', color='#555555', label="Closest point")

    ax.set_xlim(-8, 8)
    ax.set_ylim(-8, 8)
    ax.set_axisbelow(True)
    ax.legend(loc="lower left")
    ax.set_title(names[i])
    ax.grid("both", color="lightgrey", linestyle="--")
    ax.set_aspect('equal', adjustable='box')

plt.show()
