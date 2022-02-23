"""Script to test the spring components of neurites"""
import time

import numpy as np
import vedo

from neurorosettes.neurons import Neuron
from neurorosettes.simulation import Container
import neurorosettes.utilities as utilities


# Define time variables
timestep = 0.1
total_time = 40.0
pb = utilities.get_simulation_timer(total_time, timestep)

# Initialize simulation objects
container = Container(timestep=0.1,
                      viscosity=7.96,
                      grid=[-80, 80, 20])

# Populate environment with cells
container.animator.plotter.show(interactive=False, resetcam=False)

# Create a neuron with two neurites
neuron = Neuron()
neuron.create_cell(np.array([1.0, 1.0, 0.0]))
neuron.set_outgrowth_axis(np.array([1.0, 0.0, 0.0]))
neuron.create_neurites_based_on_differentiation(differentiation_grade=2)
container.register_neuron(neuron)

neuron = Neuron()
neuron.create_cell(np.array([-5.0, -15.0, 0.0]))
container.register_neuron(neuron)
neuron = Neuron()
neuron.create_cell(np.array([25.0, 6.0, 0.0]))
container.register_neuron(neuron)

container.grid.remove_cell(container.neurons[0].cell)
container.neurons[0].move_cell(np.array([-5.0, -10.0, 0.0]))
container.grid.register_cell(container.neurons[0].cell)

container.animator.plotter += vedo.Grid(sx=container.grid.representation_grid_values,
                                        sy=container.grid.representation_grid_values)

container.animator.plotter.show(interactive=False)

time.sleep(1)

# Run and plot simulation
for _ in pb.range():
    container.update_cell_positions()
    container.update_drawings()

    time.sleep(0.001)
    pb.print()

container.animator.plotter.show(interactive=True)
