"""Script to simulate the formation of a neuroblastoma rosette by cells with radial neurites (2D)"""
import numpy as np
import vedo

from neurorosettes.neurons import Neuron
from neurorosettes.simulation import Container
from neurorosettes import utilities


# Define time variables
timestep = 0.1
total_time = 10

# Initialize simulation objects
container = Container(timestep=timestep,
                      simulation_2d=True,
                      viscosity=7.96,
                      grid=[-60, 60, 20])

# Populate environment with cells
neuron1 = Neuron()
neuron1.create_cell(coordinates=np.array([2.0, 12.0, 0.0]))
container.register_neuron(neuron1)
print(neuron1.cell, container.grid.interpolate_idx(neuron1.cell.position))

neuron = Neuron()
neuron.create_cell(coordinates=np.array([25.0, 12.0, 0.0]))
container.register_neuron(neuron)

neuron = Neuron()
neuron.create_cell(coordinates=np.array([25.0, 22.0, 0.0]))
container.register_neuron(neuron)
print(neuron.cell, container.grid.get_objects_in_voxel(*container.grid.interpolate_idx(neuron.cell.position)))

neuron = Neuron()
neuron.create_cell(coordinates=np.array([-35.0, -22.0, .0]))
container.register_neuron(neuron)
print(neuron.cell, container.grid.interpolate_idx(neuron.cell.position))

print(container.grid.get_close_objects(neuron1.cell.position))
print(container.grid.get_close_objects(neuron1.cell.position))

container.animator.plotter.show(interactive=True)
