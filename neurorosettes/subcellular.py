"""This module deals with the two components of the neurons: soma cells and neurites"""
from typing import Optional, Tuple

import numpy as np
from vedo import Sphere, Spring

from neurorosettes.physics import SphereSphereInteractions, SphereCylinderInteractions, SphereMechanics, CylinderMechanics, CylinderCylinderInteractions
from neurorosettes.utilities import Animator


class CellBody:
    """Represents a cell with physical properties"""

    def __init__(self) -> None:
        self.position: np.ndarray = np.zeros(3)
        self.force_from_daughter: np.ndarray = np.zeros(3)
        self.force_from_neighbors: np.ndarray = np.zeros(3)
        self.mechanics: Optional[SphereMechanics] = None
        self.sphere: Optional[Sphere] = None

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

    def set_mechanics(self, mechanics: SphereMechanics) -> None:
        """Sets the physical properties of the cell"""
        self.mechanics = mechanics

    def set_sphere_representation(self, animator: Animator, color="red") -> None:
        """Creates and sets the sphere representation of the cell"""
        self.sphere = animator.draw_sphere(self.position, self.mechanics.radius, c=color)

    def update_representation(self) -> None:
        """Updates the sphere representation of the cell"""
        self.sphere.pos(self.position)

    def get_neighbor_force(self, neighbor: "CellBody", interaction: SphereSphereInteractions) -> np.ndarray:
        """Returns the interaction force between two cells"""
        # Compute the vector that connects the centers of the two cells
        distance_vector = np.array(self.position) - np.array(neighbor.position)
        norm = np.linalg.norm(distance_vector)
        distance_vector_normalized = distance_vector / norm

        # Calculate cell-cell adhesion forces
        adhesion_magnitude = interaction.get_adhesion_component(radius1=self.mechanics.interaction_radius,
                                                                radius2=neighbor.mechanics.interaction_radius,
                                                                distance=norm)

        adhesion_magnitude *= self.mechanics.adhesiveness * neighbor.mechanics.adhesiveness

        # Calculate cell-cell repulsion forces
        repulsion_magnitude = interaction.get_repulsion_component(radius1=self.mechanics.radius,
                                                                  radius2=neighbor.mechanics.radius,
                                                                  distance=norm)

        # Calculate resulting force
        force_magnitude = repulsion_magnitude - adhesion_magnitude

        return force_magnitude * distance_vector_normalized


class Neurite:
    """Represents a neurite with physical properties"""

    def __init__(self, proximal_point: np.ndarray, axis: np.ndarray, cylinder_mechanics: CylinderMechanics):
        """Initializes the neurite"""
        self.proximal_point: np.ndarray = proximal_point
        self.distal_point: np.ndarray = proximal_point + axis * cylinder_mechanics.default_length
        self.mechanics: CylinderMechanics = cylinder_mechanics
        self.cylinder: Optional[Tuple[Spring, Sphere]] = None
        self.force_from_daughter: Optional[np.ndarray] = None

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
        return self.mechanics.get_tension(self.current_length)

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
                                interaction: SphereCylinderInteractions) -> Tuple[np.ndarray, float]:
        """Returns the interaction force between two cells"""
        # Compute the vector that connects the centers of the two cells
        distance_vector = np.array(self.distal_point) - np.array(neighbor.position)
        norm = np.linalg.norm(distance_vector)
        distance_vector_normalized = distance_vector / norm

        # Calculate cell-cell repulsion forces
        force_magnitude, fraction = interaction.compute_interactions(cylinder_radius=self.mechanics.radius,
                                                                     cylinder_base=self.proximal_point,
                                                                     cylinder_top=self.distal_point,
                                                                     sphere_center=neighbor.position,
                                                                     sphere_radius=neighbor.mechanics.radius)

        force = force_magnitude

        return force, fraction

    def get_neurite_neighbor_force(self, neighbor: "Neurite", interaction: CylinderCylinderInteractions):
        """Returns the interaction force between two cells"""
        # Compute the vector that connects the centers of the two cells
        distance_vector = np.array(self.distal_point) - np.array(neighbor.distal_point)
        norm = np.linalg.norm(distance_vector)
        distance_vector_normalized = distance_vector / norm

        # Calculate cell-cell adhesion forces
        force, fraction = interaction.compute_interactions(self.mechanics.radius,
                                                           self.proximal_point, self.distal_point,
                                                           neighbor.mechanics.radius,
                                                           neighbor.proximal_point, neighbor.distal_point)

        return force, fraction
