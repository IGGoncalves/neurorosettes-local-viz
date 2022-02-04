"""This module deals with the neuron structure and functions"""
from typing import List

import numpy as np

from neurorosettes.physics import SphereSphereInteractions, SphereCylinderInteractions, default_potentials_repulsion, CylinderCylinderInteractions, default_potentials_adhesion
from neurorosettes.subcellular import CellBody
from neurorosettes.neurons import Neuron
from neurorosettes.utilities import Animator, get_random_unit_vector


class Container:
    """Represents the environment where neurons exist"""
    timestep: float
    sphere_interaction: SphereSphereInteractions
    sphere_cylinder_interaction: SphereCylinderInteractions
    cylinder_interaction: CylinderCylinderInteractions
    neurons: List[Neuron]
    animator: Animator

    def __init__(self, timestep: float, viscosity: float):

        self.timestep = timestep
        self.sphere_interaction = SphereSphereInteractions(default_potentials_adhesion,
                                                           default_potentials_repulsion)
        self.sphere_cylinder_interaction = SphereCylinderInteractions(default_potentials_repulsion)
        self.cylinder_interaction = CylinderCylinderInteractions(default_potentials_adhesion,
                                                                 default_potentials_repulsion)
        self.viscosity = viscosity
        self.neurons = []
        self.animator = Animator()

    def register_neuron(self, neuron: Neuron, color="red") -> None:
        """Registers a neuron and its representation into the container"""
        neuron.cell.set_sphere_representation(self.animator, color=color)
        for neurite in neuron.neurites:
            neurite.create_neurite_representation(self.animator)

        self.neurons.append(neuron)

    def update_drawings(self) -> None:
        """Updates the representations of the neurons"""
        for neuron in self.neurons:
            neuron.cell.update_representation()
            for neurite in neuron.neurites:
                neurite.update_representation()

        self.animator.plotter.show()

    def advance_cycles(self):
        for neuron in self.neurons:
            neuron.clocks.advance_clocks(self.timestep)

    def create_new_neuron(self, coordinates: np.ndarray, color="red") -> Neuron:
        """Creates a new neuron and registers it to the container"""
        new_neuron = Neuron()
        new_neuron.create_cell(coordinates)
        new_neuron.set_outgrowth_axis(get_random_unit_vector(two_dimensions=True))
        self.register_neuron(new_neuron, color=color)

        return new_neuron

    def differentiate(self) -> None:
        """Checks for neurons that are flagged for differentiation and deals with differentiation"""
        for i, neuron in enumerate(self.neurons):
            if not neuron.ready_for_differentiation or len(neuron.neurites) >= neuron.max_number_of_neurites:
                continue

            # Decide whether to create a new neurite or extend an existing one
            if neuron.neurites:
                neuron.create_secondary_neurite()
                neuron.neurites[-1].create_neurite_representation(self.animator)
            else:
                neuron.create_first_neurite()
                neuron.neurites[0].create_neurite_representation(self.animator)

            # Increase the differentiation grade and reset the differentiation signal
            neuron.clocks.differentiation_clock.increase_differentiation_grade()

    def kill(self) -> None:
        """Checks for neurons that are flagged for death and removes them from the container"""
        for neuron in self.neurons:
            if not neuron.ready_to_die:
                continue
            # Remove neuron and its representation
            self.animator.plotter -= neuron.cell.sphere
            for neurite in neuron.neurites:
                self.animator.plotter -= neurite.cylinder[0]
                self.animator.plotter -= neurite.cylinder[1]
            self.neurons.remove(neuron)

    def divide(self) -> None:
        """Checks for neurons that are flagged for division and deals with division"""
        for neuron in self.neurons:
            if not neuron.ready_for_division:
                continue
            # Create a new neuron next to the old one
            position = neuron.cell.position + get_random_unit_vector(two_dimensions=True) * neuron.cell_radius * 2.1
            color = neuron.cell.sphere.color()
            new_neuron = self.create_new_neuron(position, color)
            new_neuron.set_clocks_from_mother(neuron)
            # Update the cell cycle state of the old neuron to arrest
            neuron.clocks.cycle_clock.remove_proliferation_flag()

    @staticmethod
    def get_nearby_cells(neurons: List[Neuron]) -> List[CellBody]:
        """Returns the cells that are near to the passed cell"""
        # TODO: Build a grid to avoid long distance interactions
        return [neuron.cell for neuron in neurons]

    def get_displacement_from_force(self, force: np.ndarray) -> np.ndarray:
        """Returns a velocity from the passed force"""
        # Compute cell velocity
        velocity = force / self.viscosity
        return velocity * self.timestep

    def update_cell_positions(self) -> None:
        """Updates the cell positions and representations based on object interactions"""
        old_neurons = self.neurons.copy()
        old_cells = self.get_nearby_cells(old_neurons)

        for i, neuron in enumerate(old_neurons):
            old_neurites = neuron.neurites.copy()
            reversed_order = range(len(old_neurites) - 1, -1, -1)

            for j, old_neurite in zip(reversed_order, reversed(old_neurites)):
                force = np.zeros(3)

                # Get force from spring
                force_spring = old_neurite.get_spring_force()
                force += force_spring
                # Get force from daughter's spring
                force += old_neurite.force_from_daughter

                # Get forces from neighbors
                for k, neighbor_neuron in enumerate(old_neurons):
                    if neighbor_neuron is neuron:
                        continue
                    # Cell force
                    cell_force, fraction = old_neurite.get_cell_neighbor_force(neighbor_neuron.cell,
                                                                               self.sphere_cylinder_interaction)

                    old_neurons[k].cell.force_from_neighbors += cell_force

                    force -= cell_force * fraction
                    if j > 0:
                        neuron.neurites[j-1].force_from_daughter -= cell_force * (1 - fraction)
                    else:
                        neuron.cell.force_from_daughter -= cell_force * (1 - fraction)

                    # Neurites force
                    for m, neurite in enumerate(neighbor_neuron.neurites):
                        neurite_force, fraction = old_neurite.get_neurite_neighbor_force(neurite,
                                                                                         self.cylinder_interaction)
                        force -= neurite_force * fraction

                        if j > 0:
                            neuron.neurites[j - 1].force_from_daughter -= cell_force * (1 - fraction)
                        else:
                            neuron.cell.force_from_daughter -= cell_force * (1 - fraction)

                displacement = self.get_displacement_from_force(force)
                new_point = self.neurons[i].neurites[j].distal_point + displacement
                self.neurons[i].neurites[j].move_distal_point(new_point)

                if j > 0:
                    neuron.neurites[j-1].force_from_daughter -= force_spring
                else:
                    neuron.cell.force_from_daughter -= force_spring

                if j < len(old_neurites) - 1:
                    self.neurons[i].neurites[j+1].move_proximal_point(neuron.neurites[j].distal_point)

            force = np.zeros(3)
            for j, neighbor in enumerate(old_cells):
                if neuron.cell is neighbor:
                    continue
                force += old_cells[i].get_neighbor_force(old_cells[j], self.sphere_interaction)

            force += neuron.cell.force_from_daughter
            force += neuron.cell.force_from_neighbors

            # Convert force value to displacement to assign new position
            new_position = self.neurons[i].cell.position + self.get_displacement_from_force(force)
            self.neurons[i].cell.set_center_position(new_position)

            # Update the proximal position of the first neurite
            if neuron.neurites:
                self.neurons[i].place_neurite_on_cell_surface(neuron.neurites[0])

        for i, neuron in enumerate(old_neurons):
            old_neurons[i].cell.force_from_neighbors = np.zeros(3)
            old_neurons[i].cell.force_from_daughter = np.zeros(3)
            for neurite in neuron.neurites:
                neurite.force_from_daughter = np.zeros(3)
