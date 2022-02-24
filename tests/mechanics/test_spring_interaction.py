import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from neurorosettes import physics
from neurorosettes.subcellular import ObjectFactory

spring1 = ObjectFactory().get_neurite(proximal_position=np.array([-2.0, -5.0, 0.0]),
                                      axis=np.array([0.0, 1.0, 0.0]))

spring2 = ObjectFactory().get_neurite(proximal_position=np.array([2.0, 0.0, 0.0]),
                                      axis=np.array([0.0, 1.0, 0.0]))


point1, point2 = physics.get_cylinder_intersection(spring1.proximal_point,
                                                   spring1.distal_point,
                                                   spring2.proximal_point,
                                                   spring2.distal_point)

fig, axes = plt.subplots()
axes.add_patch(patches.Rectangle((spring1.proximal_point[0] - spring1.mechanics.radius, spring1.proximal_point[1]),
                                 width=spring1.mechanics.radius * 2,
                                 height=spring1.current_length,
                                 edgecolor="black",
                                 fill="blue"))

axes.add_patch(patches.Circle((point1[0], point1[1]), radius=spring1.mechanics.interaction_radius,
                              fill=None, alpha=1, linestyle="--"))

axes.set_xlim(-10, 10)
axes.set_ylim(-10, 10)
axes.set_axisbelow(True)
axes.grid("both", color="lightgrey")
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
