"""Script to simulate the formation of a neuroblastoma rosette by cells with radial neurites (2D)"""
import numpy as np
import vedo

from neurorosettes import physics
from neurorosettes import utilities
from neurorosettes.simulation import Container
from neurorosettes.grid import OneLevelDensityCheck

# Simulation time
timestep = 0.1
total_time = 1440

# Tissue initial configuration
tissue_size = 120.0
tissue_geometry = utilities.HexagonalTissue(tissue_size)

# Cell cycling
proliferation_rate = 0.008
death_rate = 0.0001
differentiation_rate = 0.002

# Physics
# Cell-cell interactions
sphere_sphere_adhesion = 4.0
sphere_sphere_repulsion = 40.0
smoothness_factor = 1.0
sphere_sphere_interactions = physics.get_sphere_potentials_interactions(sphere_sphere_adhesion,
                                                                        sphere_sphere_repulsion,
                                                                        smoothness_factor)
# Cell-neurite interactions
sphere_cylinder_adhesion = 0.4
sphere_cylinder_repulsion = 400.0
sphere_cylinder_interactions = physics.get_sphere_cylinder_potentials_interactions(sphere_cylinder_repulsion,
                                                                                   smoothness_factor)

# Neurite-neurite interaction
cylinder_cylinder_adhesion = 40.0
cylinder_cylinder_repulsion = 400.0
cylinder_cylinder_interactions = physics.get_cylinder_potentials_interactions(cylinder_cylinder_adhesion,
                                                                              cylinder_cylinder_repulsion,
                                                                              smoothness_factor)

# Initialize simulation objects
container = Container(timestep=timestep,
                      simulation_2d=True,
                      sphere_int=sphere_sphere_interactions,
                      sphere_cylinder_int=sphere_cylinder_interactions,
                      cylinder_int=cylinder_cylinder_interactions,
                      viscosity=7.96,
                      grid=[-160, 160, 20],
                      density_check=OneLevelDensityCheck(max_neighbors=19))

for position in tissue_geometry.get_coordinates():
    # Populate environment with cells
    neuron = container.create_new_neuron(coordinates=np.array(position))
    neuron.set_outgrowth_axis(utilities.get_random_unit_vector(two_dimensions=True))
    neuron.clocks.set_clocks(proliferation_rate, death_rate, differentiation_rate)

clock = vedo.Text2D("Simulation step: 0", pos="top right", c="black")
container.animator.plotter += clock

container.animator.plotter.show(resetcam=False, interactive=False)
pb = utilities.get_simulation_timer(total_time, timestep)

# Run and plot simulation
for t in pb.range():
    if t % 5 == 0:
        container.advance_cycles()
        container.kill()
        container.differentiate()
        container.divide()
    container.update_cell_positions()
    container.update_drawings()

    if t % 50 == 0:
        clock.text(f"Simulation step: {t}")

    pb.print()

container.animator.plotter.show(interactive=True)
