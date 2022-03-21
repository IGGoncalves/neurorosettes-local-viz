from typing import List, Tuple, Union
from dataclasses import dataclass, field
from abc import ABC, abstractmethod

import numpy as np

from neurorosettes.subcellular import CellBody, Neurite


@dataclass
class UniformGrid:
    """
    Class to store the objects inside a simulation container and their coordinates.
    
    Grids are defined as uniform objects with the same size and spacing
    for all dimensions (2D or 3D).
    """
    min: float
    """The bottom limit of the grid domain size."""
    max: float
    """The top limit of the grid domain size."""
    step: float
    """The size of the grid cells."""
    use_2d: bool = True
    """If a grid has a third dimension or is 2D."""
    grid_values: np.ndarray = field(init=False)
    """The coordinates of the grid cells."""
    idx_values: np.ndarray = field(init=False)
    """The indices of the grid cells."""
    grid: np.ndarray = field(init=False)
    """The grid, with the objects stored inside each cell."""

    def __post_init__(self) -> None:
        """Creates the empty grid object based on the initial size and spacing."""
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
        """Returns the grid cell coordinates, inluding the top limit."""
        return np.arange(self.min, self.max + 1, self.step)

    def interpolate_idx(self, position: np.ndarray) -> Tuple[int, int, int]:
        """
        Returns the indices of the cell where a given position should be stored.

        Parameters
        ----------
        position
            The position to be located inside the grid.        
        """
        idxx = np.floor(np.interp(position[0], self.grid_values, self.idx_values)).astype(int)
        idxy = np.floor(np.interp(position[1], self.grid_values, self.idx_values)).astype(int)

        if self.use_2d:
            return idxx, idxy, 0
            
        idxz = np.floor(np.interp(position[2], self.grid_values, self.idx_values)).astype(int)

        return idxx, idxy, idxz

    def register_cell(self, cell: CellBody) -> None:
        """
        Stores a cell object inside the grid and evaluates where to register it.

        Parameters
        ----------
        cell
            The cell object to be stored inside the grid.
        """
        idxx, idxy, idxz = self.interpolate_idx(cell.position)
        self.grid[idxz][idxy][idxx].append(cell)

    def register_neurite(self, neurite: Neurite) -> None:
        """
        Stores a neurite object inside the grid and evaluates where to register it.

        Parameters
        ----------
        neurite
            The neurite object to be stored inside the grid.
        """
        idxx, idxy, idxz = self.interpolate_idx(neurite.distal_point)
        self.grid[idxz][idxy][idxx].append(neurite)

    def remove_cell(self, cell: CellBody) -> None:
        """
        Removes a cell object from the grid.

        Parameters
        ----------
        cell
            The cell object to be removed from the grid.
        """
        idxx, idxy, idxz = self.interpolate_idx(cell.position)
        self.grid[idxz][idxy][idxx].remove(cell)

    def remove_neurite(self, neurite: Neurite) -> None:
        """
        Removes a neurite object from the grid.

        Parameters
        ----------
        neurite
            The neurite object to be removed from the grid.
        """
        idxx, idxy, idxz = self.interpolate_idx(neurite.distal_point)
        self.grid[idxz][idxy][idxx].remove(neurite)

    def get_objects_in_voxel(self, idxx: int, idxy: int, idxz: int  = 0) -> List[Union[CellBody, Neurite]]:
        """
        Returns the objects inside a cell of the grid.

        Parameters
        ----------
        idxx
            The index of the cell in the x component.
        idxy
            The index of the cell in the y component.
        idxz
            The index of the cell in the z component.

        Returns
        -------
        A list of the cell bodies and neurites registered in the given cell of the grid.
        """
        return self.grid[idxz][idxy][idxx]

    def get_close_objects(self, position: np.ndarray) -> List[Union[CellBody, Neurite]]:
        """
        Returns the objects inside a cell and its neighbors.

        Parameters
        ----------
        position
            The position to be evaluated.

        Returns
        -------
        A list of the cell bodies and objects inside a cell, and the
        cells that surround it in 2 or 3 dimensions.
        """
        # Find the index corresponding to the position to be evaluated
        idxx, idxy, idxz = self.interpolate_idx(position)
        neighbors = list()

        # Get the indices of the neighbors, avoiding domain boundaries
        x_neighbors = [idxx + value for value in [-1, 0, 1]
                       if 0 < idxx + value < len(self.idx_values)]

        y_neighbors = [idxy + value for value in [-1, 0, 1]
                       if 0 < idxy + value < len(self.idx_values)]

        # Disregard the z component if the grid is 2D
        if self.use_2d:
            for y in y_neighbors:
                for x in x_neighbors:
                    # Go through the cells of interest and get the objects inside
                    neighbors.extend(self.get_objects_in_voxel(x, y, 0))

            return neighbors

        z_neighbors = [idxz + value for value in [-1, 0, 1]
                       if 0 < idxz + value < len(self.idx_values)]

        # Go through the cells of interest in 3D and get the objects inside
        for z in z_neighbors:
            for y in y_neighbors:
                for x in x_neighbors:
                    neighbors.extend(self.get_objects_in_voxel(x, y, z))

        return neighbors

    def get_close_cells(self, position: np.ndarray) -> List[CellBody]:
        """
        Returns only the cell bodies inside a cell and its neighbors.

        Parameters
        ----------
        position
            The position to be evaluated.

        Returns
        -------
        A list of the cell bodies inside a cell and its neighbors.
        """
        return [neighbor for neighbor in self.get_close_objects(position)
                if isinstance(neighbor, CellBody)]

    def get_cells_in_radius(self, position: np.ndarray, radius: float) -> List[CellBody]:
        """
        Returns the objects inside a radius of interest.

        Parameters
        ----------
        position
            The position to be evaluated.
        radius
            The radius of interest.

        Returns
        -------
        A list of the cell bodies and neurites that are inside
        a radius of interest with the given position as the centre.
        """
        return [neighbor for neighbor in self.get_close_cells(position)
                if np.linalg.norm(neighbor.position - position) <= radius]


class CellDensityCheck(ABC):
    """
    Class to evaluate the density of a tissue around a neuron object.

    Multiple types of checks may be defined, such as grid-based or 
    radius-based evaluations. 
    """
    @abstractmethod
    def check_max_density(self, cell: CellBody, grid: UniformGrid) -> bool:
        """
        Returns a true boolean value if the density around a neuron surpasses the limit.
        
        Parameters
        ----------
        cell
            The cell to be evaluated.
        grid
            The Grid object where the cells and its neighbors
            are registered.
        """
        pass


class OneLevelDensityCheck(CellDensityCheck):
    """
    Class for a cell density check based on the neighboring grid cells.

    Parameters
    ----------
    max_neighbors
        The maximum numbers of neighbors inside the cell grid and
        its neighbors.
    """
    def __init__(self, max_neighbors: int = 15):
        self.max_neighbors = max_neighbors

    def check_max_density(self, cell: CellBody, grid: UniformGrid) -> bool:
        """
        Checks if the number of neighbors is higher than the maximum limit.

        Parameters
        ----------
        cell
            The cell to be evaluated.
        grid
            The Grid object where the cells and its neighbors
            are registered.
        """
        return len(grid.get_close_cells(cell.position)) >= self.max_neighbors


class TwoLevelsDensityCheck(CellDensityCheck):
    """
    Class for a cell density check based on a radius of interest.

    Starts by checking if the number of neighbors is larger than the maximum,
    as done in the OneLevelDensityCheck. If this value is surpassed, the neighbors
    inside a radius of interest are selected and quantified. It then checks if
    there are more neighbors inside this radius than the expected value.

    Parameters
    ----------
    max_outer_neighbors
        The maximum numbers of neighbors inside the cell grid and
        its neighbors.
    max_inner_neighbors
        The maximum number of neighbors inside of the radius of interest.
    radius
        The radius of interest.
    """
    def __init__(self, max_outer_neighbors: int = 15, 
                 max_inner_neighbors: int = 6, radius: float = 18.0) -> None:
        self.max_outer_neighbors = max_outer_neighbors
        self.max_inner_neighbors = max_inner_neighbors
        self.radius = radius

    def check_max_density(self, cell: CellBody, grid: UniformGrid) -> bool:
        """
        Checks if the number of neighbors is higher than the maximum limit.

        Parameters
        ----------
        cell
            The cell to be evaluated.
        grid
            The Grid object where the cells and its neighbors
            are registered.
        """
        neighbors = grid.get_close_cells(cell.position)
        if len(neighbors) >= self.max_outer_neighbors:
            return True

        cells_in_radius = grid.get_cells_in_radius(cell.position, self.radius)

        return len(cells_in_radius) >= self.max_inner_neighbors
