"""This module deals with physical interactions between objects."""
from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np


def get_sphere_overlap(distance: float, radius1: float, radius2: float) -> float:
    """Returns the overlap between two spheres"""
    # Check if cells overlap (real radius)
    overlap = radius1 + radius2 - distance
    # Return only the adhesion force if no overlap
    if overlap < 0.00001:
        return 0.0

    return overlap


@dataclass
class CylinderMechanics:
    radius: float = 0.5
    interaction_factor: float = 1.25
    spring_constant: float = 5.0
    default_length: float = 10.0

    @property
    def interaction_radius(self) -> float:
        """Returns the radius of interaction between two spheres"""
        return self.interaction_factor * self.radius

    def get_spring_tension(self, cylinder_length: float) -> float:
        """Returns the spring tension"""
        length_difference = (cylinder_length - self.default_length)
        return self.spring_constant * length_difference / self.default_length


@dataclass
class SphereMechanics:
    """Stores the physical properties that characterize a sphere"""
    radius: float
    interaction_factor: float
    adhesiveness: float

    @property
    def interaction_radius(self) -> float:
        """Returns the radius of interaction between two spheres"""
        return self.interaction_factor * self.radius


@dataclass
class Repulsion(ABC):
    """Basic class for repulsion interactions"""
    repulsion_coefficient: float

    @abstractmethod
    def compute_repulsion(self, distance: float, radius1: float, radius2: float) -> float:
        """Returns the repulsion force magnitude between two cells"""
        ...


@dataclass
class SimpleRepulsion(Repulsion):
    """Class for repulsion interactions based on CX3D"""

    def compute_repulsion(self, distance: float, radius1: float, radius2: float) -> float:
        """Returns a force proportional to the repulsion coefficient and cell overlap"""
        # Check if cells overlap (real radius)
        overlap = get_sphere_overlap(distance, radius1, radius2)

        return self.repulsion_coefficient * overlap


@dataclass
class PotentialsRepulsion(Repulsion):
    """Class for repulsion interactions based on PhysiCell"""
    smoothness_factor: float

    def compute_repulsion(self, distance: float, radius1: float, radius2: float) -> float:
        """Returns a force based on adhesion potentials"""
        # Check if cells overlap (real radius)
        overlap = get_sphere_overlap(distance, radius1, radius2)
        if overlap == 0.0:
            return 0.0
        # Compute repulsion force
        repulsion_component = (1 - distance / (radius1 + radius2))
        repulsion_component = repulsion_component ** (self.smoothness_factor + 1)

        return self.repulsion_coefficient * repulsion_component


@dataclass
class Adhesion(ABC):
    """Basic class for adhesion interactions"""
    adhesion_coefficient: float

    @abstractmethod
    def compute_adhesion(self, distance: float, radius1: float, radius2: float) -> float:
        """Returns the adhesion force magnitude between two cells"""
        ...


@dataclass
class SimpleAdhesion(Adhesion):
    """Class for adhesion interactions based on CX3D"""

    def compute_adhesion(self, distance: float, radius1: float, radius2: float) -> float:
        """Returns a force proportional to the adhesion coefficient, cell overlap and cell radii"""
        # Compute the cells' overlap, taking into account an interaction radius
        adhesion_overlap = get_sphere_overlap(distance, radius1, radius2)
        # Return no force if there is no interaction
        if adhesion_overlap == 0.0:
            return 0.0

        # Get equivalent radius
        equivalent_radius = (radius1 * radius2) / (radius1 + radius2)

        return self.adhesion_coefficient * np.sqrt(equivalent_radius * adhesion_overlap)


@dataclass
class PotentialsAdhesion(Adhesion):
    """Class for adhesion interactions based on PhysiCell"""
    smoothness_factor: int

    def compute_adhesion(self, distance: float, radius1: float, radius2: float) -> float:
        """Returns a force based on repulsion potentials"""
        # Check if cells interact (interaction radius)
        adhesion_overlap = get_sphere_overlap(distance, radius1, radius2)
        # Return no force if there is no interaction
        if adhesion_overlap == 0.0:
            return 0.0

        # Compute adhesion force
        adhesion_component = (1 - distance / (radius1 + radius2))
        adhesion_component = adhesion_component ** (self.smoothness_factor + 1)

        return self.adhesion_coefficient * adhesion_component


class SphereInteractions:
    def __init__(self, adhesion: Adhesion, repulsion: Repulsion):
        self.adhesion = adhesion
        self.repulsion = repulsion

    def get_adhesion_component(self, radius1, radius2, distance):
        return self.adhesion.compute_adhesion(distance, radius1, radius2)

    def get_repulsion_component(self, radius1, radius2, distance):
        return self.repulsion.compute_repulsion(distance, radius1, radius2)

    def get_interaction_forces(self, radius1, radius2, distance):
        force = self.adhesion.compute_adhesion(distance, radius1, radius2)
        force -= self.repulsion.compute_repulsion(distance, radius1, radius2)
        return force


class SphereCylinderInteractions:
    def __init__(self, repulsion: Repulsion) -> None:
        self.repulsion = repulsion

    @staticmethod
    def compute_intersection(sphere_center: np.ndarray, cylinder_base: np.ndarray,
                             cylinder_top: np.ndarray) -> np.ndarray:
        """Returns the intersection data (minimum distance and fraction to mother)"""
        base_to_sphere_axis = np.subtract(sphere_center, cylinder_base)
        cylinder_axis = np.subtract(cylinder_top, cylinder_base)

        dot_product = np.dot(base_to_sphere_axis, cylinder_axis)
        projection_fraction = dot_product / np.linalg.norm(cylinder_axis) ** 2

        if projection_fraction <= 0:
            return cylinder_base

        if projection_fraction >= 1:
            return cylinder_top

        return cylinder_base + cylinder_axis * projection_fraction

    def compute_interactions(self, sphere_radius: float, sphere_center: np.ndarray,
                             cylinder_radius: float, cylinder_base: np.ndarray,
                             cylinder_top: np.ndarray):

        closest_point = self.compute_intersection(sphere_center, cylinder_base, cylinder_top)
        fraction_to_mother = np.linalg.norm(np.subtract(closest_point, cylinder_base))
        fraction_to_mother /= np.linalg.norm(np.subtract(cylinder_top, cylinder_base))
        distance_vector = np.subtract(sphere_center, closest_point)
        distance = np.linalg.norm(distance_vector)

        unit_vector = distance_vector / distance
        force = self.repulsion.compute_repulsion(distance, sphere_radius, cylinder_radius)

        force *= unit_vector
        return force, fraction_to_mother


class CylinderInteractions:
    def __init__(self, adhesion: Adhesion, repulsion: Repulsion) -> None:
        self.adhesion = adhesion
        self.repulsion = repulsion

    @staticmethod
    def compute_intersection(cylinder_base_1: np.ndarray, cylinder_top_1: np.ndarray,
                             cylinder_base_2: np.ndarray, cylinder_top_2: np.ndarray):
        """Returns the intersection data (minimum distance and fraction to mother)"""
        base_axis = np.subtract(cylinder_base_1, cylinder_base_2)
        cylinder_axis_2 = np.subtract(cylinder_top_2, cylinder_base_2)
        cylinder_axis_1 = np.subtract(cylinder_top_1, cylinder_base_1)

        d1343 = np.dot(base_axis, cylinder_axis_2)
        d4321 = np.dot(cylinder_axis_2, cylinder_axis_1)
        d1321 = np.dot(base_axis, cylinder_axis_1)
        d4343 = np.dot(cylinder_axis_2, cylinder_axis_2)
        d2121 = np.dot(cylinder_axis_1, cylinder_axis_1)

        denominator = d2121 * d4343 - d4321 * d4321

        if denominator < 0.00000001:
            projection_1 = cylinder_base_1 + 0.5 * cylinder_axis_1
            projection_2 = cylinder_base_2 + 0.5 * cylinder_axis_2
            print(projection_2 - projection_1)

            return projection_1, np.subtract(projection_2, projection_1)

        numer = d1343 * d4321 - d1321 * d4343
        mua = numer / denominator
        mub = (d1343 + mua * d4321) / d4343

        if mua < 0:
            projection_1 = cylinder_base_1
        elif mua > 1:
            projection_1 = cylinder_top_1
        else:
            projection_1 = cylinder_base_1 + mua * cylinder_axis_1

        if mub < 0:
            projection_2 = cylinder_base_2
        elif mub > 1:
            projection_2 = cylinder_top_2
        else:
            projection_2 = cylinder_base_2 + mub * cylinder_axis_2

        return projection_1, np.subtract(projection_2, projection_1)

    def compute_interactions(self, cylinder_radius_1: float,
                             cylinder_base_1: np.ndarray, cylinder_top_1: np.ndarray,
                             cylinder_radius_2: float, cylinder_base_2: np.ndarray,
                             cylinder_top_2: np.ndarray):

        closest_point, distance = self.compute_intersection(cylinder_base_1, cylinder_top_1,
                                                            cylinder_base_2, cylinder_top_2)

        fraction_to_mother = np.linalg.norm(np.subtract(closest_point, cylinder_base_1))
        fraction_to_mother /= np.linalg.norm(np.subtract(cylinder_top_1, cylinder_base_1))

        distance_norm = np.linalg.norm(distance)

        unit_vector = distance / distance_norm
        force = self.repulsion.compute_repulsion(distance_norm, cylinder_radius_1, cylinder_radius_2)

        force -= self.adhesion.compute_adhesion(distance_norm, cylinder_radius_1, cylinder_radius_2)

        force *= unit_vector
        return force, fraction_to_mother


# Default interactions
def get_sphere_potentials_interactions(adhesion_coefficient, repulsion_coefficient, smoothness_factor):
    adhesion = PotentialsAdhesion(adhesion_coefficient, smoothness_factor)
    repulsion = PotentialsRepulsion(repulsion_coefficient, smoothness_factor)

    return SphereInteractions(adhesion, repulsion)


def get_sphere_cylinder_potentials_interactions(repulsion_coefficient, smoothness_factor):
    repulsion = PotentialsRepulsion(repulsion_coefficient, smoothness_factor)

    return SphereCylinderInteractions(repulsion)


def get_cylinder_potentials_interactions(adhesion_coefficient, repulsion_coefficient, smoothness_factor):
    adhesion = PotentialsAdhesion(adhesion_coefficient, smoothness_factor)
    repulsion = PotentialsRepulsion(repulsion_coefficient, smoothness_factor)

    return CylinderInteractions(adhesion, repulsion)
