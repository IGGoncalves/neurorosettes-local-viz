"""This module deals with the neuron structure and functions"""
from dataclasses import dataclass
import numpy as np

from neurorosettes.subcellular import CellBody, Neurite, ObjectFactory
from neurorosettes.clocks import CellClocks, ClocksFactory


class Neuron:
    """Represents a neuron that has a cell and neurites"""

    def __init__(self, position: np.ndarray, outgrowth_axis: np.ndarray,
                 factory: ObjectFactory, clocks: CellClocks,
                 differentiation_grade: int = 0) -> None:
        self.cell = factory.get_cell_body(position)
        self.outgrowth_axis = outgrowth_axis
        self.clocks = clocks
        self.max_number_of_neurites = 4
        self.neurites = []

        if differentiation_grade != 0:
            self.create_first_neurite(factory)
            for _ in range(1, differentiation_grade):
                self.create_secondary_neurite(factory)

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

    def create_first_neurite(self, factory: ObjectFactory) -> None:
        """Creates a neurite attached to the soma cell"""
        proximal_point = self.cell.position + self.outgrowth_axis * self.cell.mechanics.radius
        neurite = factory.get_neurite(proximal_point, self.outgrowth_axis)
        self.neurites.append(neurite)
        self.clocks.cycle_clock.cycle_block = True

    def create_secondary_neurite(self, factory: ObjectFactory) -> None:
        """Creates a neurite attached to the most recent neurite"""
        if len(self.neurites) >= self.max_number_of_neurites:
            return

        proximal_point = self.neurites[-1].distal_point
        neurite = factory.get_neurite(proximal_point, self.outgrowth_axis)
        self.neurites.append(neurite)


@dataclass
class NeuronFactory:
    objects_factory: ObjectFactory
    clocks_factory: ClocksFactory

    def create_neuron(self, coordinates, outgrowth_axis):
        clocks = self.clocks_factory.get_clocks()
        return Neuron(coordinates, outgrowth_axis, 
                      self.objects_factory, clocks)