import numpy as np
from vedo import Plotter, Sphere, Spring, ProgressBar


def get_random_position(scaling_factor: float) -> np.ndarray:
    """Returns the coordinates for a random position between -scaling_factor and scaling factor."""

    return np.array([(np.random.random()-0.5)*scaling_factor,
                     (np.random.random()-0.5)*scaling_factor,
                     (np.random.random()-0.5)*scaling_factor])


def get_random_unit_vector(two_dimensions=False) -> np.ndarray:
    """Returns a vector for a random point in a unit sphere."""
    if two_dimensions:
        vector = np.array([np.random.normal(), np.random.normal(), 0])
    else:
        vector = np.array([np.random.normal() for _ in range(3)])

    return vector / np.linalg.norm(vector)


def get_distance_unit_vector_and_norm(point1, point2):
    """Returns the direction and magnitude of a vector that connects two points."""
    distance_vector = point1 - point2
    norm = np.linalg.norm(distance_vector)
    unit_vector = distance_vector / norm

    return unit_vector, norm


def get_progress_bar(total_time: float, timestep: float) -> ProgressBar:
    """Returns a progress bar with the simulation time"""
    return ProgressBar(0, total_time / timestep, c="r")


class Animator:
    def __init__(self):
        self.plotter = Plotter(interactive=False, axes=0)
        self.set_camera()

    def set_camera(self):
        self.plotter.camera.SetPosition([0., 0., 200.])
        self.plotter.camera.SetFocalPoint([0., 0., 0.])
        self.plotter.camera.SetViewUp([0., 1., 0.])
        self.plotter.camera.SetDistance(200.)
        self.plotter.camera.SetClippingRange([180., 220.])

    def draw_spring(self, base_point: np.ndarray, top_point: np.ndarray, radius: float):
        """Plots a spring and a sphere to represent a neurite in vedo."""
        cylinder = Spring(startPoint=base_point, endPoint=top_point, r=radius)
        top_sphere = Sphere(pos=top_point, r=radius, c='b3')
        self.plotter += cylinder + top_sphere

        return cylinder, top_sphere

    def draw_sphere(self, center: np.ndarray, radius: float, **kwargs) -> Sphere:
        """Plots a sphere to represent a soma cell in vedo."""
        sphere = Sphere(pos=center, r=radius, alpha=0.6, **kwargs)
        self.plotter += sphere

        return sphere
