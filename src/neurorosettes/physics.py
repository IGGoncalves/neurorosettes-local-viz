"""This module deals with physical interactions between objects."""
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Tuple

import numpy as np


def get_distance_components(point1: np.ndarray, point2: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Returns the direction and magnitude components of a vector that connects two points.

    Parameters
    ----------
    point1: np.ndarray
        The coordinates of the first point.
    point2:np.ndarray
        The coordinates of the second point.

    Returns
    -------
    np.ndarray
        The unit vector that defines the direction between the two points.
    float
        The distance between the two points.
    """
    distance_vector = point1 - point2
    norm = np.linalg.norm(distance_vector)
    unit_vector = distance_vector / norm

    return unit_vector, norm


def get_sphere_overlap(radius1: float, radius2: float, distance: float, ) -> float:
    """
    Returns the overlap between two objects represented as spheres.

    Parameters
    ----------
    radius1: float
        The radius of the first object (represented as a sphere).
    radius2: float
        The radius of the second object (represented as a sphere).
    distance: float
        The distance between the centres of the two objects.

    Returns
    -------
    float
        The overlap between the two objects.
    """
    overlap = radius1 + radius2 - distance
    if overlap < 0.00001:
        return 0.0

    return overlap


def get_sphere_cylinder_intersection(sphere_center: np.ndarray, cylinder_base: np.ndarray,
                                     cylinder_top: np.ndarray) -> np.ndarray:
    """
    Returns the closest point on the cylinder axis to the sphere. The intersection is given by
    the dot product. Taking the dot product between the cylinder axis and the axis that connects
    the cylinder base to the sphere, the dot product gives us the projection of the cylinder-sphere
    axis on the cylinder axis. For dot products between 0 and 1, the closest point is on the
    cylinder axis. For dot products below 0 or above 1, the closest points are the base and the
    extremity of the cylinder, respectively.

    Parameters
    ----------
    sphere_center: np.ndarray
        The coordinates for the center of the sphere object.
    cylinder_base: np.ndarray
        The coordinates for the base of the cylinder object.
    cylinder_top: np.ndarray
        The coordinates for the extremity of the cylinder object.

    Returns
    -------
    np.ndarray
        The coordinates of the closest point to the sphere on the cylinder axis
    """
    base_to_sphere_axis = np.subtract(sphere_center, cylinder_base)
    cylinder_axis = np.subtract(cylinder_top, cylinder_base)

    dot_product = np.dot(base_to_sphere_axis, cylinder_axis)
    projection_fraction = dot_product / np.linalg.norm(cylinder_axis) ** 2

    if projection_fraction <= 0:
        return cylinder_base

    if projection_fraction >= 1:
        return cylinder_top

    return cylinder_base + cylinder_axis * projection_fraction


def get_cylinder_intersection(cylinder_base_1: np.ndarray,
                              cylinder_top_1: np.ndarray,
                              cylinder_base_2: np.ndarray,
                              cylinder_top_2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns the closest point between two cylinders. The cross product is used to evaluate
    if the cylinder axes are parallel. If they are, the closest points between the two objects
    are considered to be their middle points. If not, the projection of each axis on the other
    axis is computed to find the closest point on each axis.

    Parameters
    ----------
    cylinder_base_1: np.ndarray
        The coordinates for the base of the first cylinder object.
    cylinder_top_1: np.ndarray
        The coordinates for the extremity of the first cylinder object.
    cylinder_base_2: np.ndarray
        The coordinates for the base of the second cylinder object.
    cylinder_top_2: np.ndarray
        The coordinates for the extremity of the second cylinder object.

    Returns
    -------
    np.ndarray
        The closest point on the axis of the first cylinder.
    np.ndarray
        The closest point on the axis of the second cylinder.
    """
    base_axis = np.subtract(cylinder_base_1, cylinder_base_2)
    cylinder_axis_2 = np.subtract(cylinder_top_2, cylinder_base_2)
    cylinder_axis_1 = np.subtract(cylinder_top_1, cylinder_base_1)

    d1343 = np.dot(base_axis, cylinder_axis_2)
    d4321 = np.dot(cylinder_axis_2, cylinder_axis_1)
    d1321 = np.dot(base_axis, cylinder_axis_1)
    d4343 = np.dot(cylinder_axis_2, cylinder_axis_2)
    d2121 = np.dot(cylinder_axis_1, cylinder_axis_1)

    denominator = d2121 * d4343 - d4321 * d4321

    # If cylinder axes are parallel, set the closest points as the middle points
    if denominator < 0.00000001:
        projection_1 = cylinder_base_1 + 0.5 * cylinder_axis_1
        projection_2 = cylinder_base_2 + 0.5 * cylinder_axis_2

        return projection_1, projection_2

    # Calculate the projection of the cylinder axes on each other
    numer = d1343 * d4321 - d1321 * d4343
    mua = numer / denominator
    mub = (d1343 + mua * d4321) / d4343

    # Assess if the point on cylinder 1 is on the axis
    # If the projection is smaller than 0, the closest point should be the base point
    # If the projection is larger than 1, the closet point should be the extremity
    if mua < 0:
        projection_1 = cylinder_base_1
    elif mua > 1:
        projection_1 = cylinder_top_1
    else:
        projection_1 = cylinder_base_1 + mua * cylinder_axis_1

    # Assess if the point on cylinder 2 is on the axis. Use same approach as before
    #
    if mub < 0:
        projection_2 = cylinder_base_2
    elif mub > 1:
        projection_2 = cylinder_top_2
    else:
        projection_2 = cylinder_base_2 + mub * cylinder_axis_2

    return projection_1, projection_2


@dataclass
class PhysicalProperties:
    """
    Class with the mechanical properties of a physical object.

    Attributes:
    -----------
    radius: float
        The radius of the physical object.
    interaction_factor: float
        The factor used to compute physical interactions between two objects.
        Objects inside the interaction radius will be considered.
    """

    radius: float
    interaction_factor: float

    @property
    def interaction_radius(self) -> float:
        """Returns the radius of interaction of the physical object."""
        return self.interaction_factor * self.radius


@dataclass
class CylinderProperties(PhysicalProperties):
    """
    Class with the mechanical properties of a cylinder with a spring axis.

    Attributes:
    -----------
    radius: float
        The radius of the cylinder.
    interaction_factor: float
        The factor used to compute physical interactions between two spheres.
    spring_factor: float
        The spring constant used to compute the tension inside the cylinder.
    default_length: float
        The default and initial length of a cylinder.
    max_length: float
        The maximum length of a cylinder.
    """

    radius: float
    interaction_factor: float
    spring_constant: float
    default_length: float
    max_length: float

    def get_spring_tension(self, cylinder_length: float) -> float:
        """Returns the tension in the spring for a given spring length."""
        length_difference = (cylinder_length - self.default_length)
        return self.spring_constant * length_difference / self.default_length


@dataclass
class ContactForces(ABC):
    """
    Class to compute the contact forces between two objects, represented as spheres.

    Attributes
    ----------
    adhesion_coefficient: float
        The coefficient of adhesion between the two objects.
    repulsion_coefficient: float
        The coefficient of repulsion between two objects.
    """
    adhesion_coefficient: float
    repulsion_coefficient: float

    @abstractmethod
    def compute_adhesion(self, distance: float, radius1: float, radius2: float) -> float:
        """
        Returns the magnitude of the adhesion force between two objects

        Parameters
        ----------
        distance: float
            The distance between two spheres.
        radius1: float
            The radius of the first sphere object.
        radius2: float
            The radius of the second sphere object.

        Returns
        -------
        float
            The magnitude of the adhesion contact forces.
        """
        pass

    @abstractmethod
    def compute_repulsion(self, distance: float, radius1: float, radius2: float) -> float:
        """
        Returns the magnitude of the repulsion force between two objects

        Parameters
        ----------
        distance: float
            The distance between two spheres.
        radius1: float
            The radius of the first sphere object.
        radius2: float
            The radius of the second sphere object.

        Returns
        -------
        float
            The magnitude of the repulsion contact forces.
        """
        pass


@dataclass
class SimpleContact(ContactForces):
    """
    Class to compute simple contact forces between two spheres. The force components
    are proportional to the adhesion/repulsion coefficients, the overlap between the spheres and,
    in the case of adhesion forces, the equivalent radius of the spheres.
    Same approach as done in Cx3D.
    """

    def compute_adhesion(self, distance: float, radius1: float, radius2: float) -> float:
        """Returns a force proportional to the adhesion coefficient, cell overlap and cell radii"""
        # Compute the cells' overlap, taking into account an interaction radius
        adhesion_overlap = get_sphere_overlap(radius1, radius2, distance)
        # Return no force if there is no interaction
        if adhesion_overlap == 0.0:
            return 0.0

        # Get equivalent radius
        equivalent_radius = (radius1 * radius2) / (radius1 + radius2)

        return self.adhesion_coefficient * np.sqrt(equivalent_radius * adhesion_overlap)

    def compute_repulsion(self, distance: float, radius1: float, radius2: float) -> float:
        """Returns a force proportional to the repulsion coefficient and cell overlap"""
        # Check if cells overlap (real radius)
        overlap = get_sphere_overlap(radius1, radius2, distance)

        return self.repulsion_coefficient * overlap


@dataclass
class PotentialsContact(ContactForces):
    """
    Class to compute contact forces between two spheres based on potentials. The force components
    take into account the adhesion/coefficient coefficients, a smoothness factor and the distance
    between two objects.
    Same approach as done in PhysiCell.
    """
    smoothness_factor: int

    def compute_adhesion(self, distance: float, radius1: float, radius2: float) -> float:
        """Returns a force based on adhesion potentials"""
        # Check if cells interact (interaction radius)
        adhesion_overlap = get_sphere_overlap(radius1, radius2, distance)
        # Return no force if there is no interaction
        if adhesion_overlap == 0.0:
            return 0.0

        # Compute adhesion force
        adhesion_component = (1 - distance / (radius1 + radius2))
        adhesion_component = adhesion_component ** (self.smoothness_factor + 1)

        return self.adhesion_coefficient * adhesion_component

    def compute_repulsion(self, distance: float, radius1: float, radius2: float) -> float:
        """Returns a force based on repulsion potentials"""
        # Check if cells overlap (real radius)
        overlap = get_sphere_overlap(radius1, radius2, distance)
        if overlap == 0.0:
            return 0.0
        # Compute repulsion force
        repulsion_component = (1 - distance / (radius1 + radius2))
        repulsion_component = repulsion_component ** (self.smoothness_factor + 1)

        return self.repulsion_coefficient * repulsion_component


@dataclass
class ContactFactory(ABC):
    @abstractmethod
    def get_sphere_sphere_interactions(self) -> ContactForces:
        pass

    @abstractmethod
    def get_sphere_cylinder_interactions(self) -> ContactForces:
        pass

    @abstractmethod
    def get_cylinder_cylinder_interactions(self) -> ContactForces:
        pass


@dataclass
class PotentialsFactory(ContactFactory):
    sphere_sphere_adhesion: float = 4.0
    sphere_sphere_repulsion: float = 50.0
    sphere_sphere_smoothness: int = 1

    sphere_cylinder_adhesion: float = 0.0
    sphere_cylinder_repulsion: float = 100.0
    sphere_cylinder_smoothness: int = 1

    cylinder_cylinder_adhesion: float = 4.0
    cylinder_cylinder_repulsion: float = 100.0
    cylinder_cylinder_smoothness: int = 1

    def get_sphere_sphere_interactions(self) -> PotentialsContact:
        return PotentialsContact(self.sphere_sphere_adhesion,
                                 self.sphere_sphere_repulsion,
                                 self.sphere_sphere_smoothness)

    def get_sphere_cylinder_interactions(self) -> PotentialsContact:
        return PotentialsContact(self.sphere_cylinder_adhesion,
                                 self.sphere_cylinder_repulsion,
                                 self.sphere_cylinder_smoothness)

    def get_cylinder_cylinder_interactions(self) -> PotentialsContact:
        return PotentialsContact(self.cylinder_cylinder_adhesion,
                                 self.cylinder_cylinder_repulsion,
                                 self.cylinder_cylinder_smoothness)
