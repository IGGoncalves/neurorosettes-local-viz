"""This module deals with the neuron structure and functions"""
import time
from typing import List, Optional

import numpy as np

from neurorosettes.physics import SphereInteractions, SphereCylinderInteractions, CylinderInteractions
from neurorosettes.subcellular import CellBody, Neurite
from neurorosettes.neurons import Neuron
from neurorosettes.utilities import Animator, get_random_unit_vector
from neurorosettes.grid import UniformGrid, CellDensityCheck


class Container:
    """Represents the environment where neurons exist"""
    time_step: float
    simulation_2d: bool
    sphere_int: SphereInteractions
    sphere_cylinder_int: SphereCylinderInteractions
    cylinder_int: CylinderInteractions
    neurons: List[Neuron]
    animator: Animator
    grid: UniformGrid

    def __init__(self, timestep: float, simulation_2d: bool, sphere_int: SphereInteractions,
                 sphere_cylinder_int: SphereCylinderInteractions, cylinder_int: CylinderInteractions,
                 viscosity: float, grid: List[float], density_check: Optional[CellDensityCheck] = None) -> None:

        self.time_step = timestep
        self.simulation_2d = simulation_2d
        self.sphere_int = sphere_int
        self.sphere_cylinder_int = sphere_cylinder_int
        self.cylinder_int = cylinder_int
        self.viscosity = viscosity
        self.neurons = []
        self.animator = Animator()
        self.grid = UniformGrid(*grid)
        self.density_check = density_check

        if self.simulation_2d:
            self.animator.add_grid(self.grid.grid_values, self.grid.grid_values)

    def register_neuron(self, neuron: Neuron, color="red") -> None:
        """Registers a neuron and its representation into the container"""
        neuron.cell.set_sphere_representation(self.animator, color=color)
        self.grid.register_cell(neuron.cell)
        for neurite in neuron.neurites:
            neurite.create_neurite_representation(self.animator)
            self.grid.register_neurite(neurite)

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
            neuron.clocks.advance_clocks(self.time_step)

    def create_new_neuron(self, coordinates: np.ndarray, color="red") -> Neuron:
        """Creates a new neuron and registers it to the container"""
        new_neuron = Neuron()
        new_neuron.create_cell(coordinates)
        new_neuron.set_outgrowth_axis(get_random_unit_vector(two_dimensions=self.simulation_2d))

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
                neurite = neuron.neurites[-1]
                neurite.create_neurite_representation(self.animator)
                neighbors = self.grid.get_close_objects(neurite.distal_point)
                neighbor_ok = [False for _ in neighbors]

                number_of_tries = 0

                while not all(neighbor_ok) and number_of_tries < 5:
                    neighbor_ok = [False for _ in neighbors]
                    for j, neighbor in enumerate(neighbors):
                        # Check if the new neurite overlaps with another neurite
                        if isinstance(neighbor, Neurite):
                            if neighbor is neuron.neurites[-2]:
                                dot = np.dot(neurite.spring_axis/neurite.current_length,
                                             neighbor.spring_axis / neighbor.current_length)
                                if dot > 0.2:
                                    neighbor_ok[j] = True

                            else:
                                point, overlap = self.cylinder_int.compute_intersection(neurite.proximal_point,
                                                                                        neurite.distal_point,
                                                                                        neighbor.proximal_point,
                                                                                        neighbor.distal_point)

                                # If there is no overlap, mark interaction as OK
                                if np.linalg.norm(overlap) > 0.0000001:
                                    neighbor_ok[j] = True
                                else:
                                    new_point = point - neurite.spring_axis * 0.05
                                    new_length = np.linalg.norm(neurite.proximal_point - new_point)

                                    if new_length > 5.0:
                                        neurite.move_distal_point(new_point)
                                        neighbor_ok[j] = True

                        # Check if the new neurite overlaps with a cell
                        else:
                            distance = np.linalg.norm(neurite.distal_point - neighbor.position)
                            if distance >= 8:
                                neighbor_ok[j] = True
                            else:
                                unit_vector = neurite.spring_axis / np.linalg.norm(neurite.spring_axis)
                                point = neighbor.position - unit_vector * 8.5
                                if np.linalg.norm(neurite.proximal_point - point) > 5.0:
                                    neurite.move_distal_point(point)
                                    neighbor_ok[j] = True

                        if not neighbor_ok[j]:
                            new_vector = get_random_unit_vector(two_dimensions=self.simulation_2d)
                            new_vector *= neurite.mechanics.default_length
                            neurite.move_distal_point(neurite.proximal_point + new_vector)
                            break

                    number_of_tries += 1

            else:
                neuron.create_first_neurite()
                neurite = neuron.neurites[0]
                neurite.create_neurite_representation(self.animator)
                neighbors = self.grid.get_close_objects(neurite.distal_point)
                neighbor_ok = [False for _ in neighbors]

                while not all(neighbor_ok):
                    for j, neighbor in enumerate(neighbors):
                        # Check if the new neurite overlaps with another neurite
                        if isinstance(neighbor, Neurite):
                            point, overlap = self.cylinder_int.compute_intersection(neurite.proximal_point,
                                                                                    neurite.distal_point,
                                                                                    neighbor.proximal_point,
                                                                                    neighbor.distal_point)

                            # If there is no overlap, mark interaction as OK
                            if np.linalg.norm(overlap) > 0.0000001:
                                neighbor_ok[j] = True

                            else:
                                new_point = point - neurite.spring_axis * 0.05
                                new_length = np.linalg.norm(neurite.proximal_point - new_point)

                                if new_length > 5.0:
                                    neurite.move_distal_point(new_point)
                                    neighbor_ok[j] = True

                        # Check if the new neurite overlaps with a cell
                        else:
                            distance = np.linalg.norm(neurite.distal_point - neighbor.position)
                            if distance >= 8:
                                neighbor_ok[j] = True

                            else:
                                unit_vector = neurite.spring_axis / np.linalg.norm(neurite.spring_axis)
                                point = neighbor.position - unit_vector * 8.5
                                if np.linalg.norm(neurite.proximal_point - point) > 5.0:
                                    neurite.move_distal_point(point)
                                    neuron.place_neurite_on_cell_surface(neurite)
                                    neighbor_ok[j] = True

                        if not neighbor_ok[j]:
                            new_vector = get_random_unit_vector(two_dimensions=self.simulation_2d)
                            new_vector *= neurite.mechanics.default_length + neuron.cell.mechanics.radius
                            neurite.move_distal_point(neuron.cell.position + new_vector)
                            break

            if all(neighbor_ok):
                self.grid.register_neurite(neurite)
                neuron.clocks.differentiation_clock.differentiation_signal = False
            else:
                self.animator.plotter -= neurite.cylinder[0]
                self.animator.plotter -= neurite.cylinder[1]
                neuron.neurites.remove(neurite)

                neuron.clocks.differentiation_clock.differentiation_block = True

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

            if self.density_check:
                if self.density_check.check_max_density(neuron.cell, self.grid):
                    neuron.clocks.cycle_clock.remove_proliferation_flag()
                    neuron.clocks.cycle_clock.block_proliferation()
                else:
                    # Create a new neuron next to the old one
                    position = get_random_unit_vector(two_dimensions=self.simulation_2d) * neuron.cell_radius * 2.05
                    position += neuron.cell.position
                    new_neuron = self.create_new_neuron(position)
                    new_neuron.set_clocks_from_mother(neuron)
                    # Update the cell cycle state of the old neuron to arrest
                    neuron.clocks.cycle_clock.remove_proliferation_flag()
                    self.update_drawings()
            else:
                # Create a new neuron next to the old one
                position = get_random_unit_vector(two_dimensions=self.simulation_2d) * neuron.cell_radius * 2.05
                position += neuron.cell.position
                new_neuron = self.create_new_neuron(position)
                new_neuron.set_clocks_from_mother(neuron)
                # Update the cell cycle state of the old neuron to arrest
                neuron.clocks.cycle_clock.remove_proliferation_flag()
                self.update_drawings()

    def get_displacement_from_force(self, force: np.ndarray) -> np.ndarray:
        """Returns a velocity from the passed force"""
        # Compute cell velocity
        velocity = force / self.viscosity
        return velocity * self.time_step

    def update_cell_positions(self) -> None:
        """Updates the cell positions and representations based on object interactions"""
        for i, neuron in enumerate(self.neurons):
            reversed_order = range(len(neuron.neurites) - 1, -1, -1)

            for j, neurite in zip(reversed_order, reversed(neuron.neurites)):
                force = np.zeros(3)

                # Get force from spring
                force_spring = neurite.get_spring_force()
                force += force_spring
                # Get force from daughter's spring
                force += neurite.force_from_daughter
                # Get objects in the surrounding voxels
                nearby_objects = self.grid.get_close_objects(neurite.distal_point)
                nearby_cells = [nearby_object for nearby_object in nearby_objects
                                if isinstance(nearby_object, CellBody)]
                nearby_neurites = [nearby_object for nearby_object in nearby_objects
                                   if isinstance(nearby_object, Neurite)]

                # Get forces from neighbors
                for neighbor in nearby_cells:
                    if neighbor is neuron.cell:
                        continue

                    # Cell force
                    cell_force, fraction = neurite.get_cell_neighbor_force(neighbor,
                                                                           self.sphere_cylinder_int)

                    neighbor.force_from_neighbors += cell_force
                    force -= cell_force * fraction

                    # Transmit the force from cell to proximal part of the neurite
                    if j > 0:
                        neuron.neurites[j - 1].force_from_daughter -= cell_force * (1 - fraction)
                    else:
                        neuron.cell.force_from_daughter -= cell_force * (1 - fraction)

                # Neurites force
                for neighbor in nearby_neurites:
                    if neighbor in neuron.neurites:
                        continue

                    neurite_force, fraction = neurite.get_neurite_neighbor_force(neighbor,
                                                                                 self.cylinder_int)
                    force -= neurite_force * fraction

                    # Transmit the force from cell to proximal part of the neurite
                    if j > 0:
                        neuron.neurites[j - 1].force_from_daughter -= neurite_force * (1 - fraction)
                    else:
                        neuron.cell.force_from_daughter -= neurite_force * (1 - fraction)

                displacement = self.get_displacement_from_force(force)
                self.neurons[i].neurites[j].displacement = displacement

                if j > 0:
                    neuron.neurites[j - 1].force_from_daughter -= force_spring
                else:
                    neuron.cell.force_from_daughter -= force_spring

            force = np.zeros(3)
            nearby_objects = self.grid.get_close_objects(neuron.cell.position)
            nearby_cells = [nearby_object for nearby_object in nearby_objects
                            if isinstance(nearby_object, CellBody)]

            for neighbor in nearby_cells:
                if neuron.cell is neighbor:
                    continue
                force += neuron.cell.get_neighbor_force(neighbor, self.sphere_int)

            force += neuron.cell.force_from_daughter
            force += neuron.cell.force_from_neighbors

            # Convert force value to displacement to assign new position
            displacement = self.get_displacement_from_force(force)
            neuron.cell.displacement = displacement

        for neuron in self.neurons:
            neuron.cell.force_from_neighbors = np.zeros(3)
            neuron.cell.force_from_daughter = np.zeros(3)

            self.grid.remove_cell(neuron.cell)
            neuron.cell.set_center_position(neuron.cell.position + neuron.cell.displacement)
            self.grid.register_cell(neuron.cell)

            for j, neurite in enumerate(neuron.neurites):
                neurite.force_from_daughter = np.zeros(3)
                new_point = neurite.distal_point + neurite.displacement

                self.grid.remove_neurite(neurite)
                neurite.move_distal_point(new_point)
                self.grid.register_neurite(neurite)
                if j < len(neuron.neurites) - 1:
                    neuron.neurites[j + 1].move_proximal_point(neuron.neurites[j].distal_point)

            # Update the proximal position of the first neurite
            if neuron.neurites:
                neuron.place_neurite_on_cell_surface(neuron.neurites[0])
