from abc import ABC, abstractmethod

import numpy as np
from vedo import Plotter, Sphere, Spring, screenshot, Grid, Text2D


class Tissue(ABC):
    """
    Class to represent a multicellular tissue through the cell cooridnates.

    Parameters
    ----------
    use_2d
        If the cells are placed in a single plane or if
        the tissue is a 3D tissue.
    """
    def __init__(self, use_2d: bool = True):
        self.use_2d = use_2d

    @abstractmethod
    def get_coordinates(self) -> np.ndarray:
        """Returns the initial cell coordinates according to the tissue geometry"""
        pass


class RectangularTissue(Tissue):
    """
    Class to represent a tissue that is a rectangular grid.

    Parameters
    ----------
    size
        The domain bounds of the tissue.
    spacing
        The distance between cell centres.
    use_2d
        If the cells are placed in a single plane or if
        the tissue is a 3D tissue.
    """
    def __init__(self, size: float, spacing: float = 16.0, use_2d: bool = True):
        super().__init__(use_2d)
        self.size = size
        self.spacing = spacing

    def get_coordinates(self) -> np.ndarray:
        """Returns the initial cell coordinates according to the tissue geometry"""
        if self.use_2d:
            return np.array([[x, y, 0]
                             for x in np.arange(-self.size / 2, self.size / 2, self.spacing)
                             for y in np.arange(-self.size / 2, self.size / 2, self.spacing)])

        return np.array([[x, y, z]
                         for x in np.arange(-self.size / 2, self.size / 2, self.spacing)
                         for y in np.arange(-self.size / 2, self.size / 2, self.spacing)
                         for z in np.arange(-self.size / 2, self.size / 2, self.spacing)])


class HexagonalTissue(Tissue):
    """
    Class to represent a tissue that is a hexagonal grid.

    Parameters
    ----------
    size
        The domain bounds of the tissue.
    spacing
        The distance between cell centres.
    use_2d
        If the cells are placed in a single plane or if
        the tissue is a 3D tissue.
    """
    def __init__(self, size: float, spacing: float = 20.0, use_2d: bool = True):
        super().__init__(use_2d)
        self.size = size
        self.spacing = spacing

    def get_coordinates(self) -> np.ndarray:
        """Returns the initial cell coordinates according to the tissue geometry"""
        if self.use_2d:
            return np.array([[x + 12 * (i % 2), y, 0]
                             for x in np.arange(-self.size / 2, self.size / 2, self.spacing)
                             for i, y in enumerate(np.arange(-self.size / 2, self.size / 2, self.spacing))])


def get_random_position(scaling_factor: float) -> np.ndarray:
    """
    Returns the coordinates for a random position between -scaling_factor and scaling factor.
    
    Parameters
    ----------
    scaling_factor
        The range of the interval where coorinates will be
        sampled from.
    """

    return np.array([(np.random.random() - 0.5) * scaling_factor,
                     (np.random.random() - 0.5) * scaling_factor,
                     0.0])


def get_random_unit_vector(two_dimensions=False) -> np.ndarray:
    """Returns a vector for a random point in a unit sphere."""
    if two_dimensions:
        vector = np.array([np.random.normal(), np.random.normal(), 0])
    else:
        vector = np.array([np.random.normal() for _ in range(3)])

    return vector / np.linalg.norm(vector)


class Animator:
    """Class to render the simulaiton results through vedo."""
    def __init__(self):
        self.plotter = Plotter(interactive=False, axes=0, backend="ipyvtk")
        self.clock = Text2D("Simulation step: 0", pos="top right", c="black", font="Courier")
        self.plotter += self.clock

    def show(self, interactive: bool = False):
        """
        Shows the simulation results in a display window.
        
        Parameters
        ----------
        interactive
            If the plot should support user interaction.
        """
        self.plotter.show(interactive=interactive, resetcam=False)

    def update_clock(self, time_point: float) -> None:
        """
        Updates the simulation clock string with the new time point.

        Parameters
        ----------
        time_point
            The new time point to be shown on the clock string.
        """
        self.clock.text(f"Simulation step: {time_point}")

    def add_grid(self, x_grid: np.ndarray, y_grid: np.ndarray) -> None:
        """
        Adds a 2D grid object to the middle plane of the domain.

        Parameters
        ----------
        x_grid
            The x coordinates of the cells of the grid.
        y_grid
            The y coordinates of the cell of the grid.
        """
        self.plotter += Grid(sx=x_grid, sy=y_grid, c='lightgrey')

    def set_camera(self, height: float):
        """
        Places the plotter camera at a given height.

        Parameters
        ----------
        height
            The height to place the camera at.
        """
        self.plotter.camera.SetPosition([0., 0., height])
        self.plotter.camera.SetFocalPoint([0., 0., 0.])
        self.plotter.camera.SetViewUp([0., 1., 0.])

    def draw_spring(self, base_point: np.ndarray, top_point: np.ndarray, radius: float):
        """
        Plots a spring and a sphere to represent a cylinder object.

        The spring connects the two points passed as coordinates,
        and a sphere is placed centred on the top point.
        
        Parameters
        ----------
        base_point
            The coordinates for the base of the spring.
        top_point
            The coordinates for the top of the spring. 
        radius
            The radius of the cylinder object.
        """
        cylinder = Spring(startPoint=base_point, endPoint=top_point, r=radius)
        top_sphere = Sphere(pos=top_point, r=radius, c='b3')
        self.plotter += cylinder
        self.plotter += top_sphere

        return cylinder, top_sphere

    def draw_sphere(self, center: np.ndarray, radius: float, **kwargs) -> Sphere:
        """
        Plots a sphere to represent a sphere object.

        Parameters
        ----------
        center
            The coordinates for the sphere center.
        radius
            The radius of the sphere.
        """
        sphere = Sphere(pos=center, r=radius, alpha=1.0, **kwargs)
        self.plotter += sphere

        return sphere

    def save_screenshot(self, file_name: str) -> None:
        """
        Saves the current rendering to a PNG file.

        Parameters
        ----------
        file_name
            The name of the PNG file to be created.
        """
        self.plotter.screenshot(f"{file_name}.png")