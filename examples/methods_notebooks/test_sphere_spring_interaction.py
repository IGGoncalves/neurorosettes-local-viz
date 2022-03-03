import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.transforms as transforms

from neurorosettes import physics
from neurorosettes.subcellular import ObjectFactory

fig, axes = plt.subplots(1, 2, figsize=(8, 4))

names = ["Cylinder center", "Cylinder extremity"]
springs = [{"coordinates": np.array([5.5, -5.0, 0.0]), "growth": np.array([0.0, 1.0, 0]), "angle": 0},
           {"coordinates": np.array([4.0, 5.0, 0.0]), "growth": np.array([1.0, 0.0, 0]), "angle": -90}]

for i, ax in enumerate(axes):
    sphere = ObjectFactory().get_cell_body(center_position=np.array([-3.0, 0.0, 0.0]))

    spring = ObjectFactory().get_neurite(proximal_position=springs[i]["coordinates"],
                                         axis=springs[i]["growth"])

    point = physics.get_sphere_cylinder_intersection(sphere.position,
                                                     spring.proximal_point,
                                                     spring.distal_point)

    force, fraction = spring.get_cell_neighbor_force(sphere,
                                                     physics.PotentialsFactory().get_sphere_cylinder_interactions())

    ax.add_patch(patches.Circle((sphere.position[0], sphere.position[1]),
                                radius=sphere.mechanics.radius,
                                edgecolor="firebrick",
                                facecolor="salmon"))

    rect2 = patches.Rectangle((spring.proximal_point[0] - spring.mechanics.radius, spring.proximal_point[1]),
                              width=spring.mechanics.radius * 2,
                              height=spring.current_length,
                              edgecolor="darkkhaki",
                              facecolor="moccasin")

    t2 = transforms.Affine2D().rotate_deg_around(spring.proximal_point[0],
                                                 spring.proximal_point[1],
                                                 springs[i]["angle"]) + ax.transData
    rect2.set_transform(t2)
    ax.add_patch(rect2)

    ax.plot([spring.proximal_point[0], spring.distal_point[0]],
            [spring.proximal_point[1], spring.distal_point[1]],
            color="darkkhaki", linestyle="--", zorder=1)

    ax.add_patch(patches.Circle((sphere.position[0], sphere.position[1]), radius=0.2,
                                fill=True, alpha=1, facecolor="black"))

    ax.add_patch(patches.Circle((point[0], point[1]), radius=0.2,
                                fill=True, alpha=1, facecolor="black"))

    ax.add_patch(patches.Arrow(x=point[0], y=point[1], dx=force[0], dy=force[1],
                               fill=True, alpha=1, facecolor="grey"))

    ax.add_patch(patches.Arrow(x=spring.proximal_point[0], y=spring.proximal_point[1],
                               dx=force[0]*(1-fraction), dy=force[1]*(1-fraction),
                               fill=True, alpha=1, facecolor="black"))

    ax.add_patch(patches.Arrow(x=spring.distal_point[0], y=spring.distal_point[1],
                               dx=force[0] * fraction, dy=force[1] * fraction,
                               fill=True, alpha=1, facecolor="black"))

    ax.add_patch(patches.Arrow(x=sphere.position[0], y=sphere.position[1], dx=-force[0], dy=-force[1],
                               fill=True, alpha=1, facecolor="black"))

    ax.scatter([], [], s=10, marker='o', color='#555555', label="Closest point")

    ax.set_xlim(-15, 15)
    ax.set_ylim(-15, 15)
    ax.set_axisbelow(True)
    ax.legend(loc="lower left")
    ax.set_title(names[i])
    ax.grid("both", color="lightgrey", linestyle="--")
    ax.set_aspect('equal', adjustable='box')

plt.show()
