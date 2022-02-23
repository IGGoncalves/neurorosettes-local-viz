"""Script to test sphere-neurite interactions"""
import time

import numpy as np

from neurorosettes.neurons import Neuron
from neurorosettes.simulation import Container
import neurorosettes.utilities as utilities


# Define time variables
timestep = 0.1
total_time = 10.0
pb = utilities.get_simulation_timer(total_time, timestep)

# Initialize simulation objects
container = Container(timestep=0.1,
                      viscosity=7.96,
                      grid=[-80, 80, 20])

# Populate environment with cells
neuron = Neuron()
neuron.create_cell(coordinates=np.array([0.0, 0.0, 0.0]))
neuron.set_outgrowth_axis(np.array([0.0, 1.0, 0.0]))
neuron.create_neurites_based_on_differentiation(differentiation_grade=2)
container.register_neuron(neuron)

neuron = Neuron()
neuron.create_cell(coordinates=np.array([4.0, 40.0, 0.0]))
container.register_neuron(neuron)

container.animator.plotter.show(resetcam=False, interactive=False)
time.sleep(1)
# Run and plot simulation
for t in pb.range():
    container.update_cell_positions()
    container.update_drawings()
    pb.print()

container.animator.plotter.show(interactive=True)
