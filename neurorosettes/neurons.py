"""This module deals with the neuron structure and functions"""
from typing import List, Optional, Tuple

import numpy as np

from neurorosettes.physics import default_sphere_mechanics, default_cylinder_mechanics
from neurorosettes.subcellular import CellBody, Neurite
from neurorosettes.clocks import CellClocks


class Neuron:
    """Represents a neuron that has a cell and neurites"""
    cell: Optional[CellBody]
    clocks: CellClocks
    neurites: List[Neurite]
    outgrowth_axis: np.ndarray
    max_number_of_neurites: int

    def __init__(self):
        self.cell = None
        self.clocks = CellClocks()
        self.neurites = []
        self.outgrowth_axis = np.zeros(3)
        self.max_number_of_neurites = 4

    @property
    def cell_radius(self):
        return self.cell.mechanics.radius

    @property
    def ready_for_division(self):
        return self.clocks.cycle_clock.division_signal

    @property
    def ready_for_differentiation(self):
        return self.clocks.differentiation_clock.differentiation_signal

    @property
    def ready_to_die(self):
        return self.clocks.death_clock.death_signal

    def set_outgrowth_axis(self, coordinates: np.ndarray) -> None:
        """Sets the axis that neurites will follow when new neurites are created"""
        self.outgrowth_axis[0] = coordinates[0]
        self.outgrowth_axis[1] = coordinates[1]
        self.outgrowth_axis[2] = coordinates[2]

    def set_clocks_from_mother(self, mother_neuron: "Neuron"):
        clocks = mother_neuron.clocks.get_clock_rates()
        self.clocks.set_clocks(clocks[0], clocks[1], clocks[2])

    def get_neurite_position_on_cell_surface(self) -> np.ndarray:
        """Returns the coordinates on the cell surface where the first neurite should be placed"""
        connector = self.neurites[0].distal_point - self.cell.position
        unit_connector = connector / np.linalg.norm(connector)
        new_point = self.cell.position + unit_connector * self.cell.mechanics.radius

        return new_point

    def place_neurite_on_cell_surface(self, neurite: Neurite) -> None:
        """Places the neurite base on the cell surface according to the cell-distal vector"""
        neurite_attachment_coordinates = self.get_neurite_position_on_cell_surface()
        neurite.move_proximal_point(neurite_attachment_coordinates)

    def move_cell(self, coordinates: np.ndarray) -> None:
        """Moves the cell to a new position and updates the proximal point of the first neurite"""
        self.cell.set_center_position(coordinates)
        if self.neurites:
            self.place_neurite_on_cell_surface(self.neurites[0])

    def create_cell(self, coordinates: np.ndarray) -> None:
        """Creates the soma cell of the neuron at the given position"""
        self.cell = CellBody()
        self.cell.set_center_position(coordinates)
        self.cell.set_force_from_daughter(np.zeros(3))
        self.cell.set_mechanics(default_sphere_mechanics)

    def create_first_neurite(self) -> None:
        """Creates a neurite attached to the soma cell"""
        proximal_point = self.cell.position + self.outgrowth_axis * self.cell.mechanics.radius
        neurite = Neurite(proximal_point, self.outgrowth_axis, default_cylinder_mechanics)
        # TODO: This should be done automatically
        neurite.force_from_daughter = np.zeros(3)
        self.neurites.append(neurite)
        self.clocks.cycle_clock.cycle_block = True

    def create_secondary_neurite(self) -> None:
        """Creates a neurite attached to the most recent neurite"""
        if len(self.neurites) >= self.max_number_of_neurites:
            return

        proximal_point = self.neurites[-1].distal_point
        neurite = Neurite(proximal_point, self.outgrowth_axis, default_cylinder_mechanics)
        # TODO: This should be done automatically
        neurite.force_from_daughter = np.zeros(3)
        self.neurites.append(neurite)

    def create_neurites_based_on_differentiation(self, differentiation_grade: int) -> None:
        """Creates neurites based on the differentiation grade"""
        self.create_first_neurite()
        for _ in range(1, differentiation_grade):
            self.create_secondary_neurite()

