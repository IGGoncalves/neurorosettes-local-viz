"""This module deals with the neuron structure and functions"""
from pathlib import Path
from typing import List, Optional, Union
from dataclasses import  dataclass

import numpy as np
from vedo import ProgressBar

from neurorosettes.config import ConfigParser
from neurorosettes.clocks import ClocksFactory
from neurorosettes.physics import ContactFactory, PotentialsFactory
from neurorosettes.subcellular import CellBody, Neurite, ObjectFactory
from neurorosettes.neurons import Neuron, NeuronFactory
from neurorosettes.utilities import Animator, get_random_unit_vector
from neurorosettes.grid import UniformGrid, CellDensityCheck


@dataclass
class Timer:
    """Class to store the simulation time data."""
    total_time: float
    """The total time of a simulation (in minutes)."""
    step: float
    """The time between simulation points (in minutes)."""
    current_time: float = 0.0
    """The current time point of the simulation."""

    def get_progress_bar(self) -> ProgressBar:
        """Returns a progress bar with the simulation time"""
        return ProgressBar(0, self.total_time / self.step, c="r")


class Container:
    """
    Class that represents the environment where neurons exist.
    
    Parameters
    ----------
    grid
        The grid where simulation objects will be stored,
        to improve neighbor interactions.
    simulation_2d
        If the simulation is 2D or 3D.
    neuron_factory
        The factory object to be used to create new neurons.
    contact_factory
        The factory object to be used to create interactions.
    drag_coefficient
        The drag coefficient of the extracellular space.
    density_check
        Optional contact inhibition function to inhibit proliferation
        when the cell density is too high.
    """

    def __init__(self,
                 grid: UniformGrid,
                 simulation_2d: bool,
                 neuron_factory: NeuronFactory,
                 contact_factory: ContactFactory,
                 drag_coefficient: float = 10.0,
                 density_check: Optional[CellDensityCheck] = None) -> None:

        self.grid = grid
        self.simulation_2d = simulation_2d
        self.sphere_int = contact_factory.get_sphere_sphere_interactions()
        self.sphere_cylinder_int = contact_factory.get_sphere_cylinder_interactions()
        self.cylinder_int = contact_factory.get_cylinder_cylinder_interactions()
        self.neuron_factory = neuron_factory
        self.object_factory = self.neuron_factory.objects_factory
        self.drag_coefficient = drag_coefficient
        self.density_check = density_check
        self.animator = Animator()
        self.neurons = []

        if self.simulation_2d:
            self.animator.add_grid(self.grid.representation_grid_values,
                                   self.grid.representation_grid_values)

    def set_density_check(self, density_check: CellDensityCheck) -> None:
        """
        Sets the contact inhibition function to be used before proliferation.

        Parameters
        ----------
        density_check
            The contact inhibition function to be used.
        """
        self.density_check = density_check

    def register_neuron(self, neuron: Neuron, color="red") -> None:
        """
        Registers a neuron and its representation into the container.

        Parameters
        ----------
        neuron
            The new neuron to be registered.
        color
            The color of the sphere that will represent the new neuron
            in the renderings of the simulation.
        """
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

    def advance_cycles(self, time_step: float) -> None:
        """
        Updates the biological clocks of every object in the simulation.

        Parameters
        ----------
        time_step
            The time between simulation time points.
        """
        for neuron in self.neurons:
            neuron.clocks.advance_clocks(time_step)

    def create_new_neuron(self, coordinates: Union[np.ndarray, List[float]],
                          outgrowth_axis: Optional[Union[List[float], np.ndarray]] = None,
                          color="red") -> Neuron:
                          
        """
        Creates a new neuron and registers it to the container's grid.

        The new neuron is created as an undifferentiated cell body centred at
        the passed coordinates. An outgrowth axis vector can be passed to model
        neurite outgrowth along this direction.
        
        Parameters
        ----------
        coordinates
            The center position of the neuron's cell body.
        outgrowth_axis
            The direction of growth of the neuron's neurites.
        color
            The color of the new neurite in the simulation renders.
        """
        if isinstance(coordinates, list):
            coordinates = np.array(coordinates)

        if not isinstance(outgrowth_axis, np.ndarray):
            if isinstance(outgrowth_axis, list):
                outgrowth_axis = np.array(outgrowth_axis)
            else:
                outgrowth_axis = get_random_unit_vector(two_dimensions=self.simulation_2d)

        new_neuron = self.neuron_factory.create_neuron(coordinates, outgrowth_axis)
        self.register_neuron(new_neuron, color=color)

        return new_neuron

    def differentiate(self) -> None:
        """Checks for neurons that are flagged for differentiation and deals with differentiation"""
        for neuron in self.neurons:
            if not neuron.ready_for_differentiation or len(neuron.neurites) >= neuron.max_number_of_neurites:
                continue
            # Decide whether to create a new neurite or extend an existing one
            if neuron.neurites:
                neuron.create_secondary_neurite(self.object_factory)
                neurite = neuron.neurites[-1]
                neurite.create_neurite_representation(self.animator)

            else:
                neuron.create_first_neurite(self.object_factory)
                neurite = neuron.neurites[0]
                neurite.create_neurite_representation(self.animator)

            self.grid.register_neurite(neurite)
            neuron.clocks.differentiation_clock.differentiation_signal = False

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
                    neuron.clocks.cycle_clock.remove_flag()
                    neuron.clocks.cycle_clock.trigger_block()
                else:
                    # Create a new neuron next to the old one
                    position = get_random_unit_vector(two_dimensions=self.simulation_2d) * neuron.cell_radius * 2.05
                    position += neuron.cell.position
                    self.create_new_neuron(position)
                    # Update the cell cycle state of the old neuron to arrest
                    neuron.clocks.cycle_clock.remove_flag()
                    self.update_drawings()
            else:
                # Create a new neuron next to the old one
                position = get_random_unit_vector(two_dimensions=self.simulation_2d) * neuron.cell_radius * 2.05
                position += neuron.cell.position
                new_neuron = self.create_new_neuron(position)
                # Update the cell cycle state of the old neuron to arrest
                neuron.clocks.cycle_clock.remove_flag()
                self.update_drawings()

    def get_displacement_from_force(self, force: np.ndarray, time_step: float) -> np.ndarray:
        """
        Returns the displacemnt value that a force originates, based on the equation of motion.
        
        Parameters
        ----------
        force
            The force value to be converted to a displacement
        time_step
            The time passed between simulation time points.
        """
        velocity = force / self.drag_coefficient
        return velocity * time_step

    def move_cell(self, neuron: Neuron, new_coordinates: Union[np.ndarray, List[float]]) -> None:
        """
        Moves the cell to a new position and updates the proximal point of the first neurite.

        Parameters
        ----------
        neuron
            The neuron object to be moved.
        new_coordinates
            The new coordinates to be assigned to the cell body's centre.
        """
        if isinstance(new_coordinates, list):
            new_coordinates = np.array(new_coordinates)

        self.grid.remove_cell(neuron.cell)
        neuron.cell.set_center_position(new_coordinates)
        self.grid.register_cell(neuron.cell)

        if neuron.neurites:
            neuron.place_neurite_on_cell_surface(neuron.neurites[0])

    def move_neurite(self, neurite: Neurite, new_coordinates: np.ndarray) -> None:
        """
        Deals with moving a neurite's distal point and updating it on the grid.

        Parameters
        ----------
        neurite
            The neurite object to be moved.
        new_coordinates
            The new coordinates to be assigned to the neurite's distal point.
        """
        self.grid.remove_neurite(neurite)
        neurite.move_distal_point(new_coordinates)
        self.grid.register_neurite(neurite)

    def compute_displacements(self, time_step) -> None:
        """
        Computes the displacement for each object based on the resulting force.
        
        Parameters
        ----------
        time_step
            The time passed between simulation time points.
        """
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

                    neighbor.force_from_neighbors -= cell_force
                    force += cell_force * fraction

                    # Transmit the force from cell to proximal part of the neurite
                    if j > 0:
                        neuron.neurites[j - 1].force_from_daughter += cell_force * (1 - fraction)
                    else:
                        neuron.cell.force_from_daughter += cell_force * (1 - fraction)

                # Neurites force
                for neighbor in nearby_neurites:
                    if neighbor in neuron.neurites:
                        continue

                    neurite_force, fraction = neurite.get_neurite_neighbor_force(neighbor,
                                                                                 self.cylinder_int)
                    force += neurite_force * fraction

                    # Transmit the force from cell to proximal part of the neurite
                    if j > 0:
                        neuron.neurites[j - 1].force_from_daughter += neurite_force * (1 - fraction)
                    else:
                        neuron.cell.force_from_daughter += neurite_force * (1 - fraction)

                displacement = self.get_displacement_from_force(force, time_step)
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
            displacement = self.get_displacement_from_force(force, time_step)
            neuron.cell.displacement = displacement

    def update_cell_positions(self) -> None:
        """Updates the positions of all the simulation objects based on their velocity."""
        for neuron in self.neurons:
            neuron.cell.force_from_neighbors = np.zeros(3)
            neuron.cell.force_from_daughter = np.zeros(3)

            for j, neurite in enumerate(neuron.neurites):
                neurite.force_from_daughter = np.zeros(3)
                self.move_neurite(neurite, neurite.distal_point + neurite.displacement)
                if j < len(neuron.neurites) - 1:
                    neuron.neurites[j + 1].move_proximal_point(neuron.neurites[j].distal_point)

            # Update the proximal position of the first neurite
            self.move_cell(neuron, neuron.cell.position + neuron.cell.displacement)

    def solve_mechanics(self, time_step) -> None:
        """
        Solves the mechanical interactions and updates the neurons' positions.

        Goes through each object and computes the resulting force acting on
        it, then gets the object's velocity based on the equation of motion.
        When all of the objects are checked, the positions are updated based
        on the calculated velocity.

        Parameters
        ----------
        time_step
            The time passed between simulation time points.
        """
        self.compute_displacements(time_step)
        self.update_cell_positions()


class Simulation:
    """
    Class to create and run a simulation.

    Parameters
    ----------
    timer
        The structure to store the time data of the simulation.
    container
        The structure to store the spatial data of the simulation.
    """
    def __init__(self, timer: Timer, container: Container):
        self.timer = timer
        self.container = container

    def run(self) -> None:
        """Runs the entire simulation by solving the mechanics at each time point."""
        sim_time = self.timer.get_progress_bar()

        for t in sim_time.range():
            self.container.advance_cycles(self.timer.step)
            self.container.kill()
            self.container.differentiate()
            self.container.divide()
            # Solve interactions and draw the new object positions
            self.container.solve_mechanics(self.timer.step)
            self.container.update_drawings()

            # Update the simulation time on the simulation window
            if t % 10 == 0:
                self.container.animator.update_clock(t)

            # Print time to the console as a progressbar
            self.timer.current_time += self.timer.step
            sim_time.print()

    @classmethod
    def from_file(cls, config_path: Union[Path, str]) -> "Simulation":
        """
        Initializes a Simulation object from a YAML config file.

        Parameters
        ----------
        config_path
            The path to the YAML file config file.
        """
        if not isinstance(config_path, Path):
            config_path = Path(config_path)

        parser = ConfigParser(config_path)

        timer = Timer(**parser.get_time_data())
        grid = UniformGrid(**parser.get_domain_data())
        status_2d = parser.get_2d_status()
        drag = parser.get_drag_coefficient()

        number_of_neurites = parser.get_max_number_of_neurites()
        objects = ObjectFactory(**parser.get_objects_data())
        clocks = ClocksFactory(**parser.get_clocks_data())

        interactions_data = parser.get_interactions_data()
        interactions_type = interactions_data.pop("type")
        if interactions_type == "potentials":
            interactions = PotentialsFactory(**interactions_data)
        else:
            interactions = PotentialsFactory(**interactions_data)

        container = Container(grid=grid,
                              simulation_2d=status_2d,
                              neuron_factory=NeuronFactory(number_of_neurites, objects, clocks),
                              contact_factory=interactions,
                              drag_coefficient=drag)

        return Simulation(timer, container)
