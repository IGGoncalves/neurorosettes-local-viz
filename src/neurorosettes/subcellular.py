"""This module deals with the two components of the neurons: soma cells and neurites"""
import time
from typing import Optional, Tuple
from dataclasses import dataclass

import numpy as np
from vedo import Sphere, Spring

from neurorosettes import physics
from neurorosettes.utilities import Animator


class CellBody:
    """Represents a cell with physical properties"""

    def __init__(self, position: np.ndarray, mechanics: physics.PhysicalProperties) -> None:
        self.position = position
        self.mechanics = mechanics
        self.force_from_daughter = np.zeros(3)
        self.force_from_neighbors = np.zeros(3)
        self.displacement = np.zeros(3)
        self.sphere = None

    def set_center_position(self, coordinates: np.ndarray) -> None:
        """Sets the cell position"""
        self.position[0] = coordinates[0]
        self.position[1] = coordinates[1]
        self.position[2] = coordinates[2]

    def set_force_from_daughter(self, force: np.ndarray) -> None:
        """Sets the force transmitted by daughter"""
        self.force_from_daughter[0] = force[0]
        self.force_from_daughter[1] = force[1]
        self.force_from_daughter[2] = force[2]

    def set_mechanics(self, mechanics: physics.PhysicalProperties) -> None:
        """Sets the physical properties of the cell"""
        self.mechanics = mechanics

    def set_sphere_representation(self, animator: Animator, color="red") -> None:
        """Creates and sets the sphere representation of the cell"""
        self.sphere = animator.draw_sphere(self.position, self.mechanics.radius, c=color)

    def update_representation(self) -> None:
        """Updates the sphere representation of the cell"""
        self.sphere.pos(self.position)

    def get_neighbor_force(self, neighbor: "CellBody", interaction: physics.ContactForces) -> np.ndarray:
        """Returns the interaction force between two cells"""
        # Compute the vector that connects the centers of the two cells
        distance_vector, norm = physics.get_distance_components(neighbor.position, self.position)

        # Calculate cell-cell adhesion forces
        magnitude = interaction.compute_adhesion(distance=norm, radius1=self.mechanics.interaction_radius,
                                                 radius2=neighbor.mechanics.interaction_radius)

        magnitude -= interaction.compute_repulsion(norm, self.mechanics.radius, neighbor.mechanics.radius)

        return magnitude * distance_vector


class Neurite:
    """Represents a neurite with physical properties"""

    def __init__(self, proximal_point: np.ndarray, axis: np.ndarray, cylinder_mechanics: physics.CylinderProperties):
        """Initializes the neurite"""
        self.proximal_point: np.ndarray = proximal_point
        self.distal_point: np.ndarray = proximal_point + axis * cylinder_mechanics.default_length
        self.mechanics = cylinder_mechanics
        self.cylinder: Optional[Tuple[Spring, Sphere]] = None
        self.force_from_daughter = np.zeros(3)
        self.displacement = np.zeros(3)

    def set_force_from_daughter(self, force):
        """Sets the force transmitted by daughter"""
        self.force_from_daughter = force

    def create_neurite_representation(self, animator: Animator):
        """Creates a spring+sphere representation of the neurite"""
        spring, sphere = animator.draw_spring(self.proximal_point, self.distal_point, self.mechanics.radius)
        self.cylinder = (spring, sphere)

    def update_representation(self):
        """Updates the representation of the neurite"""
        self.cylinder[0].stretch(self.proximal_point, self.distal_point)
        self.cylinder[1].pos(self.distal_point)

    @property
    def tension(self) -> float:
        """Returns the tension of the spring"""
        return self.mechanics.get_spring_tension(self.current_length)

    @property
    def spring_axis(self) -> np.ndarray:
        """Returns the vector that defines the spring"""
        return self.distal_point - self.proximal_point

    @property
    def current_length(self) -> float:
        """Returns the current length of the spring"""
        return np.linalg.norm(self.spring_axis)

    def move_distal_point(self, coordinates: np.ndarray) -> None:
        """Moves the distal point of the spring"""
        self.distal_point = coordinates

    def move_proximal_point(self, coordinates) -> None:
        """Moves the proximal point of the spring"""
        self.proximal_point = coordinates

    def get_growth_force(self, magnitude: float):
        """Returns the force created by neurite growth"""
        return magnitude * self.spring_axis

    def get_spring_force(self) -> np.ndarray:
        """Returns the force created by spring tension"""
        return -self.tension / self.current_length * self.spring_axis

    def get_cell_neighbor_force(self, neighbor: CellBody,
                                interaction: physics.ContactForces) -> Tuple[np.ndarray, float]:
        """Returns the interaction force between two cells"""
        point = physics.get_sphere_cylinder_intersection(sphere_center=neighbor.position,
                                                         cylinder_base=self.proximal_point,
                                                         cylinder_top=self.distal_point)

        distance_to_point = np.linalg.norm(np.subtract(point, self.proximal_point))
        fraction_to_mother = distance_to_point / self.current_length

        distance_vector, norm = physics.get_distance_components(neighbor.position, point)

        # Calculate cell-cell adhesion forces
        magnitude = interaction.compute_adhesion(distance=norm, radius1=self.mechanics.interaction_radius,
                                                 radius2=neighbor.mechanics.interaction_radius)

        magnitude -= interaction.compute_repulsion(norm, self.mechanics.radius, neighbor.mechanics.radius)

        return magnitude * distance_vector, fraction_to_mother

    def get_neurite_neighbor_force(self, neighbor: "Neurite", interaction: physics.ContactForces):
        """Returns the interaction force between two cells"""
        point1, point2 = physics.get_cylinder_intersection(cylinder_base_1=self.proximal_point,
                                                           cylinder_top_1=self.distal_point,
                                                           cylinder_base_2=neighbor.proximal_point,
                                                           cylinder_top_2=neighbor.distal_point)

        distance_to_point = np.linalg.norm(np.subtract(point1, self.proximal_point))
        fraction_to_mother = distance_to_point / self.current_length

        distance_vector, norm = physics.get_distance_components(point2, point1)

        # Calculate cell-cell adhesion forces
        magnitude = interaction.compute_adhesion(distance=norm, radius1=self.mechanics.interaction_radius,
                                                 radius2=neighbor.mechanics.interaction_radius)

        magnitude -= interaction.compute_repulsion(distance=norm, radius1=self.mechanics.radius,
                                                   radius2=neighbor.mechanics.radius)

        return magnitude * distance_vector, fraction_to_mother


@dataclass
class ObjectFactory:
    cell_radius: float = 8.0
    cell_interaction_factor: float = 1.25
    neurite_radius: float = 0.5
    neurite_interaction_factor: float = 2.5
    neurite_spring_constant: float = 10.0
    neurite_default_length: float = 10.0
    neurite_max_length: float = 15.0

    def get_cell_body(self, center_position: np.ndarray) -> CellBody:
        cell_mechanics = physics.PhysicalProperties(self.cell_radius,
                                                    self.cell_interaction_factor)

        return CellBody(center_position, cell_mechanics)

    def get_neurite(self, proximal_position: np.ndarray, axis: np.ndarray) -> Neurite:
        cylinder_mechanics = physics.CylinderProperties(self.neurite_radius,
                                                        self.neurite_interaction_factor,
                                                        self.neurite_spring_constant,
                                                        self.neurite_default_length,
                                                        self.neurite_max_length)

        return Neurite(proximal_position, axis, cylinder_mechanics)
