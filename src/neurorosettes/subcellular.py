"""This module deals with the two components of the neurons: soma cells and neurites"""
from typing import Tuple
from dataclasses import dataclass

import numpy as np

from neurorosettes import physics
from neurorosettes.utilities import Animator


class CellBody:
    """
    Class to represent a cell body with physical properties.

    Parameters
    ----------
    position
        The coordinates of the cell body centre.
    mechanics
        The physical properties of the sphere that defines
        the cell body object.
    """

    def __init__(
        self, position: np.ndarray, mechanics: physics.PhysicalProperties
    ) -> None:
        self.position = position
        self.mechanics = mechanics
        self.force = np.zeros(3)
        self.force_from_daughter = np.zeros(3)
        self.force_from_neighbors = np.zeros(3)
        self.displacement = np.zeros(3)
        self.sphere = None

    def set_center_position(self, coordinates: np.ndarray) -> None:
        """
        Sets the cell position based on the given coordinates.

        Parameters
        ----------
        coordinates
            The cocordinates of the cell body centre.
        """
        self.position[0] = coordinates[0]
        self.position[1] = coordinates[1]
        self.position[2] = coordinates[2]

    def set_force_from_daughter(self, force: np.ndarray) -> None:
        """
        Sets the force transmitted by the cell's daughter, if it exists.

        Parameters
        ----------
        force
            The force components of the force transmitted by the
            daughter object.
        """
        self.force_from_daughter[0] = force[0]
        self.force_from_daughter[1] = force[1]
        self.force_from_daughter[2] = force[2]

    def set_mechanics(self, mechanics: physics.PhysicalProperties) -> None:
        """Sets the physical properties of the cell"""
        self.mechanics = mechanics

    def set_sphere_representation(self, animator: Animator, color="red") -> None:
        """Creates and sets the sphere representation of the cell"""
        self.sphere = animator.draw_sphere(
            self.position, self.mechanics.radius, c=color
        )

    def update_representation(self) -> None:
        """Updates the sphere representation of the cell"""
        self.sphere.pos(self.position)

    def get_neighbor_force(
        self, neighbor: "CellBody", interaction: physics.ContactForces
    ) -> np.ndarray:
        """Returns the interaction force between two cells"""
        # Compute the vector that connects the centers of the two cells
        distance_vector, norm = physics.get_distance_components(
            neighbor.position, self.position
        )

        # Calculate cell-cell adhesion forces
        magnitude = interaction.compute_adhesion(
            distance=norm,
            radius1=self.mechanics.interaction_radius,
            radius2=neighbor.mechanics.interaction_radius,
        )

        magnitude -= interaction.compute_repulsion(
            norm,
            self.mechanics.radius,
            neighbor.mechanics.radius,
        )

        return magnitude * distance_vector


class Neurite:
    """
    Class to represent a neurite with physical properties.

    Parameters
    ----------
    proximal_point
        The coordinates of the neurite's proximal point.
    axis
        The axis of growth used to calculate the neurite's distal point.
    cylinder_mechanics
        The physical properties of the sphere that defines
        the cell body object.
    """

    def __init__(
        self,
        proximal_point: np.ndarray,
        axis: np.ndarray,
        cylinder_mechanics: physics.CylinderProperties,
    ):
        """Initializes the neurite"""
        self.proximal_point = proximal_point
        self.distal_point = proximal_point + axis * cylinder_mechanics.default_length
        self.mechanics = cylinder_mechanics
        self.cylinder = None
        self.force = np.zeros(3)
        self.force_from_daughter = np.zeros(3)
        self.displacement = np.zeros(3)

    def set_force_from_daughter(self, force):
        """
        Sets the force transmitted by the cell's daughter, if it exists.

        Parameters
        ----------
        force
            The force components of the force transmitted by the
            daughter object.
        """
        self.force_from_daughter = force

    def create_neurite_representation(self, animator: Animator):
        """Creates a spring+sphere representation of the neurite."""
        spring, sphere = animator.draw_spring(
            self.proximal_point, self.distal_point, self.mechanics.radius
        )
        self.cylinder = (spring, sphere)

    def update_representation(self):
        """Updates the Animator representation of the neurite."""
        self.cylinder[0].stretch(self.proximal_point, self.distal_point)
        self.cylinder[1].pos(self.distal_point)

    @property
    def tension(self) -> float:
        """Returns the tension of the spring."""
        return self.mechanics.get_spring_tension(self.current_length)

    @property
    def spring_axis(self) -> np.ndarray:
        """Returns the vector that defines the spring."""
        return self.distal_point - self.proximal_point

    @property
    def current_length(self) -> float:
        """Returns the current length of the spring."""
        return np.linalg.norm(self.spring_axis)

    def move_distal_point(self, coordinates: np.ndarray) -> None:
        """Moves the distal point of the spring to the given position."""
        self.distal_point = coordinates

    def move_proximal_point(self, coordinates) -> None:
        """Moves the proximal point of the spring to the given position."""
        self.proximal_point = coordinates

    def get_growth_force(self, magnitude: float):
        """Returns the force created by neurite growth"""
        return magnitude * self.spring_axis

    def get_spring_force(self) -> np.ndarray:
        """Returns the force created by spring tension"""
        return -self.tension / self.current_length * self.spring_axis

    def get_cell_neighbor_force(
        self, neighbor: CellBody, interaction: physics.ContactForces
    ) -> Tuple[np.ndarray, float]:
        """
        Returns the interaction force between the neurite and a CellBody.

        Parameters
        ----------
        neighbor
            The CellBody with which the object is interacting.
        interactions
            The contact force to be used to calculate the interaction force.

        Returns
        -------
        np.ndarray
            The components of the interaction force between the two objects.
        float
            The fraction of the force to be applied to the distal point.
        """
        # Get the point on the cylinder axis closest to the sphere
        point = physics.get_sphere_cylinder_intersection(
            center=neighbor.position, base=self.proximal_point, top=self.distal_point
        )

        # Calculate the distance between the two closest points
        distance_to_point = np.linalg.norm(np.subtract(point, self.proximal_point))
        fraction = distance_to_point / self.current_length
        distance_vector, norm = physics.get_distance_components(
            neighbor.position, point
        )

        # Calculate cell-cell adhesion forces
        magnitude = interaction.compute_adhesion(
            distance=norm,
            radius1=self.mechanics.interaction_radius,
            radius2=neighbor.mechanics.interaction_radius,
        )

        magnitude -= interaction.compute_repulsion(
            norm,
            self.mechanics.radius,
            neighbor.mechanics.radius,
        )

        return magnitude * distance_vector, fraction

    def get_neurite_neighbor_force(
        self, neighbor: "Neurite", interaction: physics.ContactForces
    ) -> Tuple[np.ndarray, float]:
        """
        Returns the interaction force between the neurite and another Neurite.

        Parameters
        ----------
        neighbor
            The other Neurite with which the object is interacting.
        interactions
            The contact force to be used to calculate the interaction force.

        Returns
        -------
        np.ndarray
            The components of the interaction force between the two objects.
        float
            The fraction of the force to be applied to the distal point.
        """
        # Get the closes point on the cylinder axes to the other cylinder
        point1, point2 = physics.get_cylinder_intersection(
            base_1=self.proximal_point,
            top_1=self.distal_point,
            base_2=neighbor.proximal_point,
            top_2=neighbor.distal_point,
        )

        # print("Internal check", point1, point2)

        # Calculate the distance between the two closest points
        distance_to_point = np.linalg.norm(np.subtract(self.proximal_point, point1))
        fraction = distance_to_point / self.current_length
        distance_vector, norm = physics.get_distance_components(point2, point1)

        # Calculate cell-cell adhesion forces
        magnitude = interaction.compute_adhesion(
            distance=norm,
            radius1=self.mechanics.interaction_radius,
            radius2=neighbor.mechanics.interaction_radius,
        )

        magnitude -= interaction.compute_repulsion(
            distance=norm,
            radius1=self.mechanics.radius,
            radius2=neighbor.mechanics.radius,
        )

        return magnitude * distance_vector, fraction


@dataclass
class ObjectFactory:
    """
    Helper class to create instances of cell bodies and neurties based on defined properties.
    """

    cell_radius: float
    """The cell radius."""
    cell_interaction_factor: float
    """The factor used to calculate the cell's interaction radius."""
    neurite_radius: float
    """The neurite radius."""
    neurite_interaction_factor: float
    """The factor used to calculate the neurite's interaction radius."""
    neurite_spring_constant: float
    """The neurite's spring constant."""
    neurite_default_length: float
    """The initial default length of a new neuritte."""

    def get_cell_body(self, center_position: np.ndarray) -> CellBody:
        """
        Returns a CellBody object centred on the given position.

        Parameters
        ----------
        center_position
            The coordinates of the cell body's centre.

        Returns
        -------
        A CellBody with the given coordinates and the chosen properties.
        """
        cell_mechanics = physics.PhysicalProperties(
            self.cell_radius, self.cell_interaction_factor
        )

        return CellBody(center_position, cell_mechanics)

    def get_neurite(
        self,
        proximal_position: np.ndarray,
        axis: np.ndarray,
        restriction_factor: float = 0.8,
    ) -> Neurite:
        """
        Returns a Neurite object placed on the given position, pointing towards the axis.

        Parameters
        ----------
        proximal_position
            The coordinates of the neurite's proximal point.
        axis
            The axis of growth of the new neurite.
        restriction_factor
            The restriction factor used to compute the direction of a new neurite
            (the neurite will follow the axis with a given degree of freedom, defined
            by the restriction factor). Should be between 0 and 1.

        Returns
        -------
        A Neurite with the given coordinates and the chosen properties.
        """
        if np.linalg.norm(axis) > 1.0:
            axis = physics.normalize_vector(axis)

        if not (0 <= restriction_factor <= 1):
            raise ValueError(
                "The neurite axis restriction factor should be between 0 and 1."
            )

        axis_angle = np.arctan2(axis[1], axis[0])
        temp_angle = (
            np.pi
            * (
                (1 - restriction_factor) * np.random.uniform()
                + (-0.5 + restriction_factor / 2)
            )
            + axis_angle
        )
        new_axis = np.asarray([np.cos(temp_angle), np.sin(temp_angle), 0])

        cylinder_mechanics = physics.CylinderProperties(
            self.neurite_radius,
            self.neurite_interaction_factor,
            self.neurite_spring_constant,
            self.neurite_default_length,
        )

        return Neurite(proximal_position, new_axis, cylinder_mechanics)
