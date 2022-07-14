"""This module deals with the neuron structure and functions"""
from dataclasses import dataclass
import numpy as np

from neurorosettes.subcellular import Neurite, ObjectFactory
from neurorosettes.clocks import CellClocks, ClocksFactory
from neurorosettes.physics import normalize_vector


class Neuron:
    """
    Class to represent a neuron as a combination of a cell body and neurites.

    Parameters
    ----------
    position
        The center position of the neuron's cell body.
    outgrowth_axis
        The direction of growth of the neuron's neurites.
    factory
        The factory object to be used to create new neuron components.
    clocks
        The rates of the neuron's biological processes.
    max_number_of_neurites
        The total number of neurites a neuron can have before
        proliferation is blocked.
    differentiation_grade
        The starting differentiation stage of the neurite, represented
        by the initial number of neurites.
    """

    def __init__(
        self,
        position: np.ndarray,
        outgrowth_axis: np.ndarray,
        factory: ObjectFactory,
        clocks: CellClocks,
        max_number_of_neurites: int = 4,
        differentiation_grade: int = 0,
    ) -> None:
        self.cell = factory.get_cell_body(position)
        self.outgrowth_axis = outgrowth_axis
        self.clocks = clocks
        self.max_number_of_neurites = max_number_of_neurites
        self.neurites = []

        if differentiation_grade != 0:
            self.create_first_neurite(factory)
            for _ in range(1, differentiation_grade):
                self.create_secondary_neurite(factory)

    @property
    def cell_radius(self):
        """The radius of the neuron's cell body."""
        return self.cell.mechanics.radius

    @property
    def ready_for_division(self):
        """If the neuron is ready to proliferate."""
        return self.clocks.cycle_clock.signal

    @property
    def ready_for_differentiation(self):
        """If the neuron is ready to differentiate."""
        return self.clocks.differentiation_clock.signal

    @property
    def ready_to_die(self):
        """If the neuron is ready to die."""
        return self.clocks.death_clock.signal

    def set_outgrowth_axis(self, coordinates: np.ndarray) -> None:
        """
        Sets the axis that neurites will follow when new neurites are created.

        Prameters
        ----------
        outgrowth_axis
            The direction of growth of the neuron's neurites.
        """
        if np.linalg.norm(coordinates) > 1.0:
            coordinates = normalize_vector(coordinates)
            
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
        """
        Places the neurite base on the cell surface according to the cell-distal vector.

        Parameters
        ----------
        neurite
            The neurite to be placed on the cell body's surface.
        """
        neurite_attachment_coordinates = self.get_neurite_position_on_cell_surface()
        neurite.move_proximal_point(neurite_attachment_coordinates)

    def create_first_neurite(self, factory: ObjectFactory) -> None:
        """
        Creates a neurite attached to the soma cell

        Parameters
        ----------
        factory
            The factory object to be used to create new neuron components.
        """
        proximal_point = (
            self.cell.position + self.outgrowth_axis * self.cell.mechanics.radius
        )
        neurite = factory.get_neurite(proximal_point, self.outgrowth_axis)
        self.neurites.append(neurite)
        self.clocks.cycle_clock.cycle_block = True

    def create_secondary_neurite(self, factory: ObjectFactory) -> None:
        """
        Creates a neurite attached to the most recent neurite.

        Parameters
        ----------
        factory
            The factory object to be used to create new neuron components.
        """
        if len(self.neurites) >= self.max_number_of_neurites:
            return

        proximal_point = self.neurites[-1].distal_point
        neurite = factory.get_neurite(proximal_point, self.outgrowth_axis)
        self.neurites.append(neurite)


@dataclass
class NeuronFactory:
    """
    Helper class to create instances of neurons based on defined properties.
    """

    max_number_of_neurites: int
    """The maximum number of neurites that a neuron can have."""
    objects_factory: ObjectFactory
    """The ObjectFactory used to create new object instances."""
    clocks_factory: ClocksFactory
    """The ClocksFactory used to create new clock instances."""

    def create_neuron(
        self, coordinates: np.ndarray, outgrowth_axis: np.ndarray
    ) -> Neuron:
        """Returns a Neuron object placed in the given coordinates."""
        clocks = self.clocks_factory.get_clocks()
        return Neuron(
            coordinates,
            outgrowth_axis,
            self.objects_factory,
            clocks,
            self.max_number_of_neurites,
        )
