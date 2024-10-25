"""This module deals with physical interactions between objects."""
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Tuple

import numpy as np


def get_distance_components(
    point1: np.ndarray, point2: np.ndarray
) -> Tuple[np.ndarray, float]:
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
    if norm < 0.0000001:
        norm = 0.0000001

    unit_vector = distance_vector / norm

    return unit_vector, norm


def normalize_vector(vector: np.ndarray) -> np.ndarray:
    """
    Returns the input array normalized to a unit vector.

    Parameters
    ----------
    vector
        The array to be normalized.

    Returns
    -------
    np.ndarray
        The normalized vector.
    """
    return vector / np.linalg.norm(vector)


def get_sphere_overlap(radius1: float, radius2: float, distance: float) -> float:
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


def get_sphere_cylinder_intersection(
    center: np.ndarray, base: np.ndarray, top: np.ndarray
) -> np.ndarray:
    """
    Returns the closest point on the cylinder axis to the sphere.

    The intersection is given by
    the dot product. Taking the dot product between the cylinder axis and the axis that connects
    the cylinder base to the sphere, the dot product gives us the projection of the cylinder-sphere
    axis on the cylinder axis. For dot products between 0 and 1, the closest point is on the
    cylinder axis. For dot products below 0 or above 1, the closest points are the base and the
    extremity of the cylinder, respectively.

    Parameters
    ----------
    center: np.ndarray
        The coordinates for the center of the sphere object.
    base: np.ndarray
        The coordinates for the base of the cylinder object.
    top: np.ndarray
        The coordinates for the extremity of the cylinder object.

    Returns
    -------
    np.ndarray
        The coordinates of the closest point to the sphere on the cylinder axis
    """
    base_to_sphere_axis = np.subtract(center, base)
    cylinder_axis = np.subtract(top, base)

    dot_product = np.dot(base_to_sphere_axis, cylinder_axis)
    projection_fraction = dot_product / np.linalg.norm(cylinder_axis) ** 2

    if projection_fraction <= 0:
        return base

    if projection_fraction >= 1:
        return top

    return base + cylinder_axis * projection_fraction


def get_cylinder_intersection(
    base_1: np.ndarray, top_1: np.ndarray, base_2: np.ndarray, top_2: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns the closest point between two cylinders.

    The cross product is used to evaluate
    if the cylinder axes are parallel. If they are, the closest points between the two objects
    are considered to be their middle points. If not, the projection of each axis on the other
    axis is computed to find the closest point on each axis.

    Parameters
    ----------
    base_1: np.ndarray
        The coordinates for the base of the first cylinder object.
    top_1: np.ndarray
        The coordinates for the extremity of the first cylinder object.
    base_2: np.ndarray
        The coordinates for the base of the second cylinder object.
    top_2: np.ndarray
        The coordinates for the extremity of the second cylinder object.

    Returns
    -------
    np.ndarray
        The closest point on the axis of the first cylinder.
    np.ndarray
        The closest point on the axis of the second cylinder.
    """
    cylinder_axis_2 = np.subtract(top_2, base_2)
    cylinder_axis_1 = np.subtract(top_1, base_1)

    unit_vector_2, norm2 = get_distance_components(top_2, base_2)
    unit_vector_1, norm1 = get_distance_components(top_1, base_1)

    cross = np.cross(unit_vector_1, unit_vector_2)

    denominator = np.linalg.norm(cross) ** 2

    # If cylinder axes are parallel, set the closest points as the middle points
    if denominator < 0.000001:
        d0 = np.dot(unit_vector_1, (base_2 - base_1))
        d1 = np.dot(unit_vector_1, (top_2 - base_1))

        if d0 <= 0 >= d1:
            if np.absolute(d0) < np.absolute(d1):
                return base_1, base_2
            return base_1, top_2

        elif d0 >= norm1 <= d1:
            if np.absolute(d0) < np.absolute(d1):
                return top_1, base_2
            return top_1, top_2

        p1 = base_1 + 0.5 * cylinder_axis_1
        p2 = base_2 + 0.5 * cylinder_axis_2

        return p1, p2

    t = base_2 - base_1
    detA = np.linalg.det([t, unit_vector_2, cross])
    detB = np.linalg.det([t, unit_vector_1, cross])

    t0 = detA / denominator
    t1 = detB / denominator

    # Closest point on segment 1
    pA = base_1 + (unit_vector_1 * t0)
    # Closest point on segment B 2
    pB = base_2 + (unit_vector_2 * t1)

    # Clamp projections
    if t0 < 0:
        pA = base_1
    elif t0 > norm1:
        pA = top_1

    if t1 < 0:
        pB = base_2
    elif t1 > norm2:
        pB = top_2

    # Clamp projection A
    if (t0 < 0) or (t0 > norm1):
        dot = np.dot(unit_vector_2, (pA - base_2))
        if dot < 0:
            dot = 0
        elif dot > norm2:
            dot = norm2
        pB = base_2 + (unit_vector_2 * dot)

    # Clamp projection B
    if (t1 < 0) or (t1 > norm2):
        dot = np.dot(unit_vector_1, (pB - base_1))
        if dot < 0:
            dot = 0
        elif dot > norm1:
            dot = norm1
        pA = base_1 + (unit_vector_1 * dot)

    return pA, pB


class PhysicalProperties:
    """
    Class with the mechanical properties of a physical object.

    Parameters
    ----------
    radius
        The radius of the physical object.
    interaction_factor
        The factor used to calculate the radius of interaction.
    """

    def __init__(self, radius: float, interaction_factor: float) -> None:
        self.radius = radius
        self.interaction_factor = interaction_factor

    @property
    def radius(self):
        """
        Returns the radius of the object.
        Raises an error if a negative or non-numerical value is set.
        """
        return self._radius

    @radius.setter
    def radius(self, radius: float):
        if not isinstance(radius, (int, float)):
            raise TypeError("Object radius should be a numerical value.")
        if radius < 0.0:
            raise ValueError("Object radius should be positive.")
        self._radius = radius

    @property
    def interaction_factor(self):
        """
        Returns the interaction factor of the object.
        Raises an error if a negative or non-numerical value is set.
        """
        return self._interaction_factor

    @interaction_factor.setter
    def interaction_factor(self, interaction_factor: float):
        if not isinstance(interaction_factor, (int, float)):
            raise TypeError("Interaction factor should be a numerical value.")
        if interaction_factor < 0.0:
            raise ValueError("Interaction factor should be positive.")
        self._interaction_factor = interaction_factor

    @property
    def interaction_radius(self) -> float:
        """Returns the radius of interaction of the physical object."""
        return self.interaction_factor * self.radius


class SphereProperties(PhysicalProperties):
    """
    Class with the mechanical properties of a sphere.

    Parameters
    ----------
    radius
        The radius of the physical object.
    interaction_factor
        The factor used to calculate the radius of interaction.
    """

    def __init__(self, radius: float, interaction_factor: float):
        super().__init__(radius, interaction_factor)


class CylinderProperties(PhysicalProperties):
    """
    Class with the mechanical properties of a cylinder with a spring axis.

    Parameters
    ----------
    radius
        The radius of the physical object.
    interaction_factor
        The factor used to calculate the radius of interaction.
    """

    def __init__(
        self,
        radius: float,
        interaction_factor: float,
        spring_constant: float,
        default_length: float,
    ) -> None:
        super().__init__(radius, interaction_factor)
        self.spring_constant = spring_constant
        self.default_length = default_length

    def get_spring_tension(self, cylinder_length: float) -> float:
        """Returns the tension in the spring for a given spring length."""
        length_difference = cylinder_length - self.default_length
        return self.spring_constant * length_difference / self.default_length


@dataclass
class ContactForces(ABC):
    """Class to compute the contact forces between two objects, represented as spheres."""

    adhesion_coefficient: float
    repulsion_coefficient: float

    @abstractmethod
    def compute_adhesion(
        self, distance: float, radius1: float, radius2: float
    ) -> float:
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
    def compute_repulsion(
        self, distance: float, radius1: float, radius2: float
    ) -> float:
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

    def compute_adhesion(
        self, distance: float, radius1: float, radius2: float
    ) -> float:
        """Returns a force proportional to the adhesion coefficient, cell overlap and cell radii"""
        # Compute the cells' overlap, taking into account an interaction radius
        adhesion_overlap = get_sphere_overlap(radius1, radius2, distance)
        # Return no force if there is no interaction
        if adhesion_overlap == 0.0:
            return 0.0

        # Get equivalent radius
        equivalent_radius = (radius1 * radius2) / (radius1 + radius2)

        return self.adhesion_coefficient * np.sqrt(equivalent_radius * adhesion_overlap)

    def compute_repulsion(
        self, distance: float, radius1: float, radius2: float
    ) -> float:
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

    def compute_adhesion(
        self, distance: float, radius1: float, radius2: float
    ) -> float:
        """Returns a force based on adhesion potentials"""
        # Check if cells interact (interaction radius)
        adhesion_overlap = get_sphere_overlap(radius1, radius2, distance)
        # Return no force if there is no interaction
        if adhesion_overlap == 0.0:
            return 0.0

        # Compute adhesion force
        adhesion_component = 1 - distance / (radius1 + radius2)
        adhesion_component = adhesion_component ** (self.smoothness_factor + 1)

        return self.adhesion_coefficient * adhesion_component

    def compute_repulsion(
        self, distance: float, radius1: float, radius2: float
    ) -> float:
        """Returns a force based on repulsion potentials"""
        # Check if cells overlap (real radius)
        overlap = get_sphere_overlap(radius1, radius2, distance)
        if overlap == 0.0:
            return 0.0
        # Compute repulsion force
        repulsion_component = 1 - distance / (radius1 + radius2)
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
class SimpleFactory(ContactFactory):
    sphere_sphere_adhesion: float
    sphere_sphere_repulsion: float

    sphere_cylinder_adhesion: float
    sphere_cylinder_repulsion: float

    cylinder_cylinder_adhesion: float
    cylinder_cylinder_repulsion: float

    def get_sphere_sphere_interactions(self) -> ContactForces:
        return SimpleContact(self.sphere_sphere_adhesion, self.sphere_sphere_repulsion)

    def get_sphere_cylinder_interactions(self) -> ContactForces:
        return SimpleContact(
            self.sphere_cylinder_adhesion, self.sphere_cylinder_repulsion
        )

    def get_cylinder_cylinder_interactions(self) -> ContactForces:
        return SimpleContact(
            self.cylinder_cylinder_adhesion, self.cylinder_cylinder_repulsion
        )


@dataclass
class PotentialsFactory(ContactFactory):
    sphere_sphere_adhesion: float
    sphere_sphere_repulsion: float
    sphere_sphere_smoothness: int

    sphere_cylinder_adhesion: float
    sphere_cylinder_repulsion: float
    sphere_cylinder_smoothness: int

    cylinder_cylinder_adhesion: float
    cylinder_cylinder_repulsion: float
    cylinder_cylinder_smoothness: int

    def get_sphere_sphere_interactions(self) -> PotentialsContact:
        return PotentialsContact(
            self.sphere_sphere_adhesion,
            self.sphere_sphere_repulsion,
            self.sphere_sphere_smoothness,
        )

    def get_sphere_cylinder_interactions(self) -> PotentialsContact:
        return PotentialsContact(
            self.sphere_cylinder_adhesion,
            self.sphere_cylinder_repulsion,
            self.sphere_cylinder_smoothness,
        )

    def get_cylinder_cylinder_interactions(self) -> PotentialsContact:
        return PotentialsContact(
            self.cylinder_cylinder_adhesion,
            self.cylinder_cylinder_repulsion,
            self.cylinder_cylinder_smoothness,
        )
