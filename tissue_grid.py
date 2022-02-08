"""Script to simulate the formation of a neuroblastoma rosette by cells with radial neurites (2D)"""
import glob
import os
import time
from PIL import Image

import numpy as np
import vedo

from neurorosettes.neurons import Neuron
from neurorosettes.simulation import Container
from neurorosettes.physics import default_sphere_interactions
from neurorosettes import utilities


def create_tissue() -> list:
    """Returns a list of randomly distributed positions for a given number of cells"""
    coordinates = [np.array([x + 8 * (i % 2), y, 0])
                   for x in range(-60, 61, 17)
                   for i, y in enumerate(range(-45, 46, 15))]
    return coordinates


# Define time variables
timestep = 0.1
total_time = 400
pb = utilities.get_progress_bar(total_time, timestep)

# Initialize simulation objects
container = Container(timestep=timestep,
                      viscosity=7.96,
                      grid=[-160, 160, 20])

for position in create_tissue():
    # Populate environment with cells
    neuron = Neuron()
    neuron.create_cell(coordinates=np.array(position))
    neuron.set_outgrowth_axis(utilities.get_random_unit_vector(two_dimensions=True))
    neuron.clocks.set_clocks(0.005, 0.0001, 0.01)
    container.register_neuron(neuron)

container.animator.plotter.show(resetcam=False, interactive=False)
time.sleep(5)

# Run and plot simulation
for i, t in enumerate(pb.range()):
    container.advance_cycles()
    container.kill()
    container.differentiate()
    container.divide()
    container.update_cell_positions()
    container.update_drawings()
    #if i % 1000:
    #    vedo.io.screenshot(f"{str(i).zfill(5)}.png")
    pb.print()

vedo.io.screenshot(f"differentiation_on_7_proliferation_on.png")
container.animator.plotter.close()
