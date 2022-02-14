"""Script to simulate the formation of a neuroblastoma rosette by cells with radial neurites (2D)"""
import numpy as np
import vedo

from neurorosettes.neurons import Neuron
from neurorosettes.simulation import Container
from neurorosettes import utilities


def create_tissue() -> list:
    """Returns a list of randomly distributed positions for a given number of cells"""
    coordinates = [utilities.get_random_position(100) for _ in range(50)]
    return coordinates


# Define time variables
timestep = 0.1
total_time = 1440
pb = utilities.get_progress_bar(total_time, timestep)

clock = vedo.Text2D("Simulation step: 0", pos="top right", c="black")

# Initialize simulation objects
container = Container(timestep=timestep,
                      simulation_2d=False,
                      viscosity=7.96,
                      grid=[-160, 160, 20])

for position in create_tissue():
    # Populate environment with cells
    neuron = Neuron()
    neuron.create_cell(coordinates=np.array(position))
    neuron.set_outgrowth_axis(utilities.get_random_unit_vector(two_dimensions=False))
    neuron.clocks.set_clocks(0.00008, 0.0001, 0.002)
    container.register_neuron(neuron)

container.animator.plotter += clock
container.animator.plotter.show(resetcam=False, interactive=False)

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
