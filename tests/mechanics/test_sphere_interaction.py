import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.transforms as transforms

from neurorosettes import physics
from neurorosettes.subcellular import ObjectFactory

fig, ax = plt.subplots()

sphere1 = ObjectFactory().get_cell_body(center_position=np.array([-6.5, -4.5, 0.0]))

sphere2 = ObjectFactory().get_cell_body(center_position=np.array([6.5, 8.5, 0.0]))

ax.add_patch(patches.Circle((sphere1.position[0], sphere1.position[1]),
                            radius=sphere1.mechanics.radius,
                            edgecolor="firebrick",
                            facecolor="salmon"))

ax.add_patch(patches.Circle((sphere2.position[0], sphere2.position[1]),
                            radius=sphere1.mechanics.radius,
                            edgecolor="firebrick",
                            facecolor="salmon"))

ax.add_patch(patches.Circle((sphere1.position[0], sphere1.position[1]),
                            radius=sphere1.mechanics.interaction_radius,
                            fill=None, alpha=1, facecolor="black", linestyle="--"))

ax.add_patch(patches.Circle((sphere2.position[0], sphere2.position[1]),
                            radius=sphere2.mechanics.interaction_radius,
                            fill=None, alpha=1, facecolor="black", linestyle="--"))

ax.add_patch(patches.Circle((sphere1.position[0], sphere1.position[1]),
                            radius=0.2,
                            fill=True, alpha=1, facecolor="black"))

ax.add_patch(patches.Circle((sphere2.position[0], sphere2.position[1]),
                            radius=0.2,
                            fill=True, alpha=1, facecolor="black"))

ax.plot([], [], color='#555555', linestyle="--", label="Interaction radius")

ax.set_xlim(-20, 20)
ax.set_ylim(-20, 20)
ax.set_axisbelow(True)
ax.legend(loc="lower left")
ax.grid("both", color="lightgrey", linestyle="--")
ax.set_aspect('equal', adjustable='box')

plt.show()
