from abc import ABC, abstractmethod

import numpy as np
from vedo import Plotter, Sphere, Spring, ProgressBar, Grid, Text2D


class Tissue(ABC):
    def __init__(self, use_2d: bool = True):
        self.use_2d = use_2d

    @abstractmethod
    def get_coordinates(self) -> np.ndarray:
        """Returns the initial cell coordinates according to the tissue geometry"""
        pass


class RectangularTissue(Tissue):
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
    """Returns the coordinates for a random position between -scaling_factor and scaling factor."""

    return np.array([(np.random.random() - 0.5) * scaling_factor,
                     (np.random.random() - 0.5) * scaling_factor,
                     (np.random.random() - 0.5) * scaling_factor])


def get_random_unit_vector(two_dimensions=False) -> np.ndarray:
    """Returns a vector for a random point in a unit sphere."""
    if two_dimensions:
        vector = np.array([np.random.normal(), np.random.normal(), 0])
    else:
        vector = np.array([np.random.normal() for _ in range(3)])

    return vector / np.linalg.norm(vector)


class Animator:
    def __init__(self):
        self.plotter = Plotter(interactive=False, axes=0)
        self.clock = Text2D("Simulation step: 0", pos="top right", c="black", font="Courier")
        self.plotter += self.clock

    def show(self, interactive: bool = False):
        self.plotter.show(interactive=interactive, resetcam=False)

    def update_clock(self, time_point: float) -> None:
        self.clock.text(f"Simulation step: {time_point}")

    def add_grid(self, x_grid, y_grid):
        self.plotter += Grid(sx=x_grid, sy=y_grid, c='lightgrey')

    def set_camera(self, height: float):
        self.plotter.camera.SetPosition([0., 0., height])
        self.plotter.camera.SetFocalPoint([0., 0., 0.])
        self.plotter.camera.SetViewUp([0., 1., 0.])

    def draw_spring(self, base_point: np.ndarray, top_point: np.ndarray, radius: float):
        """Plots a spring and a sphere to represent a neurite in vedo."""
        cylinder = Spring(startPoint=base_point, endPoint=top_point, r=radius)
        top_sphere = Sphere(pos=top_point, r=radius, c='b3')
        self.plotter += cylinder
        self.plotter += top_sphere

        return cylinder, top_sphere

    def draw_sphere(self, center: np.ndarray, radius: float, **kwargs) -> Sphere:
        """Plots a sphere to represent a soma cell in vedo."""
        sphere = Sphere(pos=center, r=radius, alpha=1.0, **kwargs)
        self.plotter += sphere

        return sphere
