from typing import List, Tuple, Union
from dataclasses import dataclass, field
from abc import ABC, abstractmethod

import numpy as np

from neurorosettes.subcellular import CellBody, Neurite


@dataclass
class UniformGrid:
    min: float
    max: float
    step: float
    use_2d: bool = True
    grid_values: np.ndarray = field(init=False)
    idx_values: np.ndarray = field(init=False)
    grid: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        self.grid_values = np.arange(self.min, self.max, self.step)
        self.idx_values = np.arange(self.grid_values.shape[0])

        z_shape = 1 if self.use_2d else self.grid_values.shape[0]

        self.grid = np.empty(shape=[z_shape,
                                    self.grid_values.shape[0],
                                    self.grid_values.shape[0]], dtype=object)

        self.grid[...] = [[[[] for _ in range(self.grid.shape[2])]
                           for _ in range(self.grid.shape[1])]
                          for _ in range(self.grid.shape[0])]

    @property
    def representation_grid_values(self) -> np.ndarray:
        return np.arange(self.min, self.max + 1, self.step)

    def interpolate_idx(self, position: np.ndarray) -> Tuple[int, int, int]:
        idxx = np.floor(np.interp(position[0], self.grid_values, self.idx_values)).astype(int)
        idxy = np.floor(np.interp(position[1], self.grid_values, self.idx_values)).astype(int)

        if self.use_2d:
            return idxx, idxy, 0
            
        idxz = np.floor(np.interp(position[2], self.grid_values, self.idx_values)).astype(int)

        return idxx, idxy, idxz

    def register_cell(self, cell: CellBody) -> None:
        idxx, idxy, idxz = self.interpolate_idx(cell.position)
        self.grid[idxz][idxy][idxx].append(cell)

    def register_neurite(self, neurite: Neurite) -> None:
        idxx, idxy, idxz = self.interpolate_idx(neurite.distal_point)
        self.grid[idxz][idxy][idxx].append(neurite)

    def remove_cell(self, cell: CellBody) -> None:
        idxx, idxy, idxz = self.interpolate_idx(cell.position)
        self.grid[idxz][idxy][idxx].remove(cell)

    def remove_neurite(self, neurite: Neurite) -> None:
        idxx, idxy, idxz = self.interpolate_idx(neurite.distal_point)
        self.grid[idxz][idxy][idxx].remove(neurite)

    def get_objects_in_voxel(self, idxx: int, idxy: int, idxz: int  = 0) -> List[Union[CellBody, Neurite]]:
        return self.grid[idxz][idxy][idxx]

    def get_close_objects(self, position: np.ndarray) -> List[Union[CellBody, Neurite]]:
        idxx, idxy, idxz = self.interpolate_idx(position)
        neighbors = list()

        x_neighbors = [idxx + value for value in [-1, 0, 1]
                       if 0 < idxx + value < len(self.idx_values)]

        y_neighbors = [idxy + value for value in [-1, 0, 1]
                       if 0 < idxy + value < len(self.idx_values)]

        if self.use_2d:
            for y in y_neighbors:
                for x in x_neighbors:
                    neighbors.extend(self.get_objects_in_voxel(x, y, 0))

            return neighbors

        z_neighbors = [idxz + value for value in [-1, 0, 1]
                       if 0 < idxz + value < len(self.idx_values)]

        for z in z_neighbors:
            for y in y_neighbors:
                for x in x_neighbors:
                    neighbors.extend(self.get_objects_in_voxel(x, y, z))

        return neighbors

    def get_close_cells(self, position: np.ndarray) -> List[CellBody]:
        return [neighbor for neighbor in self.get_close_objects(position)
                if isinstance(neighbor, CellBody)]

    def get_cells_in_radius(self, position: np.ndarray, radius: float) -> List[CellBody]:
        return [neighbor for neighbor in self.get_close_cells(position)
                if np.linalg.norm(neighbor.position - position) <= radius]


class CellDensityCheck(ABC):
    @abstractmethod
    def check_max_density(self, cell: CellBody, grid: UniformGrid) -> bool:
        pass


class OneLevelDensityCheck(CellDensityCheck):
    def __init__(self, max_neighbors: int = 15):
        self.max_neighbors = max_neighbors

    def check_max_density(self, cell: CellBody, grid: UniformGrid) -> bool:
        return len(grid.get_close_cells(cell.position)) >= self.max_neighbors


class TwoLevelsDensityCheck(CellDensityCheck):
    def __init__(self, max_outer_neighbors: int = 15, max_inner_neighbors: int = 6, radius: float = 18.0) -> None:
        self.max_outer_neighbors = max_outer_neighbors
        self.max_inner_neighbors = max_inner_neighbors
        self.radius = radius

    def check_max_density(self, cell: CellBody, grid: UniformGrid) -> bool:
        neighbors = grid.get_close_cells(cell.position)
        if len(neighbors) >= self.max_outer_neighbors:
            return True

        cells_in_radius = grid.get_cells_in_radius(cell.position, self.radius)

        return len(cells_in_radius) >= self.max_inner_neighbors
