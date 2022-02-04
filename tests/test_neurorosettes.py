"""Script to test the spring components of neurites"""
import time

import numpy as np

from neurorosettes.neurons import Neuron
from neurorosettes.simulation import Container
from neurorosettes.physics import default_sphere_interactions
import neurorosettes.utilities as utilities

# Define time variables
timestep = 0.1
total_time = 10.0
pb = utilities.get_progress_bar(total_time, timestep)

# Initialize simulation objects
container = Container(timestep=0.1,
                      viscosity=7.96)

# Populate environment with cells
container.animator.plotter.show(interactive=False, resetcam=False)

# Create a neuron with two neurites
neuron = Neuron()
neuron.create_cell(np.zeros(3))
neuron.set_outgrowth_axis(np.array([1.0, 0.0, 0.0]))
neuron.create_neurites_based_on_differentiation(differentiation_grade=2)
container.register_neuron(neuron)

# Plot the current state of the simulation
container.animator.plotter.show(interactive=False)
time.sleep(1)

# Move the cell to test if the springs pull the parts together
container.neurons[0].move_cell(np.array([0.0, -20.0, 0.0]))

# Plot the current state of the simulation
container.update_drawings()
container.animator.plotter.show(interactive=False, resetcam=False)
time.sleep(1)

# Run and plot simulation
for _ in pb.range():
    container.advance_cycles()
    container.kill()
    container.divide()
    container.differentiate()
    container.update_cell_positions()
    container.update_drawings()
    pb.print()

container.animator.plotter.show(interactive=True)
