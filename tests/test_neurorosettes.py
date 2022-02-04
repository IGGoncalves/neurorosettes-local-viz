"""Script to test the spring components of neurites"""
from dataclasses import dataclass, field
from typing import List, Union

import numpy as np
import vedo

from neurorosettes.subcellular import CellBody, Neurite
from neurorosettes.neurons import Neuron
from neurorosettes.simulation import Container
import neurorosettes.utilities as utilities


@dataclass
class UniformGrid:
    min: float
    max: float
    step: float
    grid_values: np.ndarray = field(init=False)
    idx_values: np.ndarray = field(init=False)
    grid: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        self.grid_values = np.arange(self.min, self.max, self.step)
        self.idx_values = np.arange(self.grid_values.shape[0])
        self.grid = np.empty(shape=[self.grid_values.shape[0], self.grid_values.shape[0]], dtype=object)
        self.grid[...] = [[[] for _ in range(self.grid.shape[0])] for _ in range(self.grid.shape[1])]

    @property
    def representation_grid_valuesd(self) -> np.ndarray:
        return np.arange(self.min, self.max+1, self.step)

    def interpolate_idx(self, position: np.ndarray) -> List[int]:
        idxx = np.floor(np.interp(position[0], self.grid_values, self.idx_values)).astype(int)
        idxy = np.floor(np.interp(position[1], self.grid_values, self.idx_values)).astype(int)

        return idxx, idxy

    def register_cell(self, cell: CellBody) -> None:
        idxx, idxy = self.interpolate_idx(cell.position)
        self.grid[idxy][idxx].append(cell)

    def register_neurite(self, neurite: Neurite) -> None:
        idxx, idxy = self.interpolate_idx(neurite.distal_point)
        self.grid[idxy][idxx].append(neurite)

    def remove_cell(self, cell: CellBody) -> None:
        idxx, idxy = self.interpolate_idx(cell.position)
        self.grid[idxy][idxx].remove(cell)

    def remove_neurite(self, neurite: Neurite) -> None:
        idxx, idxy = self.interpolate_idx(neurite.distal_point)
        self.grid[idxy][idxx].remove(neurite)

    def get_objects_in_voxel(self, idxx, idxy) -> List[Union[CellBody, Neurite]]:
        return self.grid[idxy][idxx]

    def get_close_objects(self, position: np.ndarray) -> List[Union[CellBody, Neurite]]:
        idxx, idxy = self.interpolate_idx(position)
        neighbors = self.get_objects_in_voxel(idxx, idxy)

        if idxx-1 > self.idx_values[0]:
            neighbors.extend(self.get_objects_in_voxel(idxx-1, idxy))
            if idxy > self.idx_values[0]:
                neighbors.extend(self.get_objects_in_voxel(idxx-1, idxy-1))
            if idxy < self.idx_values[-1]+1:
                neighbors.extend(self.get_objects_in_voxel(idxx-1, idxy+1))

        if idxx+1 < self.idx_values[-1]+1:
            neighbors.extend(self.get_objects_in_voxel(idxx+1, idxy))
            if idxy > self.idx_values[0]:
                neighbors.extend(self.get_objects_in_voxel(idxx+1, idxy-1))
            if idxy < self.idx_values[-1]+1:
                neighbors.extend(self.get_objects_in_voxel(idxx+1, idxy+1))

        return neighbors


# Define time variables
timestep = 0.1
total_time = 10.0
pb = utilities.get_progress_bar(total_time, timestep)

# Initialize simulation objects
container = Container(timestep=0.1,
                      viscosity=7.96)

domain = UniformGrid(-40, 40, 20)

# Populate environment with cells
container.animator.plotter.show(interactive=False, resetcam=False)

# Create a neuron with two neurites
neuron = Neuron()
neuron.create_cell(np.array([1.0, 1.0, 0.0]))
neuron.set_outgrowth_axis(np.array([1.0, 0.0, 0.0]))
neuron.create_neurites_based_on_differentiation(differentiation_grade=2)
container.register_neuron(neuron)

# Create a neuron with two neurites
neuron2 = Neuron()
neuron2.create_cell(np.array([-15.0, -18.0, 0.0]))
container.register_neuron(neuron2)

domain.register_cell(neuron.cell)
domain.register_cell(neuron2.cell)
domain.register_neurite(neuron.neurites[0])
domain.register_neurite(neuron.neurites[1])

print(domain.get_close_objects(neuron.neurites[0].distal_point))

container.animator.plotter += vedo.Grid(sx=domain.representation_grid_valuesd, 
                                        sy=domain.representation_grid_valuesd)

container.animator.plotter.show(interactive=True)
