"""Script to perform a parameter study on the cell mechanics parameters that regulate the tissue architecture."""
from typing import Dict
from pathlib import Path
import json
from itertools import count

from neurorosettes.simulation import Simulation, Container
from neurorosettes.utilities import HexagonalTissue
from neurorosettes.grid import OneLevelDensityCheck

# Cell-cell interactions
proliferation_rate = [0.0000000000000001, 0.0007]
differentiation_rate = [0.0006, 0.0000000000000001]
number_of_neurites = [10, 3]

# Replicates for each condition
NUMBER_OF_REPLICATES = 2

# Contact inhibition function
DENSITY_CHECK = OneLevelDensityCheck(max_neighbors=19)


def set_density_check(container: Container) -> None:
    """Adds the contact inhibition function to the simulation."""
    container.set_density_check(DENSITY_CHECK)


def create_tissue(container: Container) -> None:
    """Creates and registers new neurons in the simulation world."""
    # Initial cell positions
    tissue = HexagonalTissue(size=220, spacing=20).get_coordinates()
    # Create a neuron with two neurites
    for cell in tissue:
        container.create_new_neuron(coordinates=cell)


def simulation(sim_world: Simulation) -> None:
    """Adds any additional functions to the container and runs the simulation."""
    # Add additional functions to the simulation (i.g., contact inhibition)
    set_density_check(sim_world.container)
    # Create initial configuration
    create_tissue(sim_world.container)
    # Plot the current state of the simulation
    sim_world.container.animator.set_camera(height=450.0)
    sim_world.container.animator.show()
    # Run the simulation to check if springs work
    sim_world.run()


def create_params_dict(prol_rate: float, diff_rate: float, max_neurites: int) -> Dict[str, float]:
    """Returns a dictionary with the name of the varied parameters and corresponding values."""
    return {
        "proliferation_rate": prol_rate,
        "differentiation_rate": diff_rate,
        "maximum_number_of_neurites": max_neurites
    }

def main():
    # Simple counter to make sure that output files are stored in unique folders
    sim_counter = count(start=0)

    # Loop through the interaction combinations
    for max_neurites, prol_rate in zip(number_of_neurites,proliferation_rate):
        for diff_rate in differentiation_rate:
            # Store the current iteration's parameters in a JSON file
            # First create a dict with this information
            parameters = create_params_dict(prol_rate, diff_rate, max_neurites)

            # Then store the data in a JSON file
            folder_name = next(sim_counter)
            filepath = Path(f"output/{folder_name}/params.json")
            filepath.parent.mkdir(parents=True, exist_ok=True)

            with open(filepath, "w") as outfile:
                json.dump(parameters, outfile)

            for i in range(NUMBER_OF_REPLICATES):
                # Load the data from the configuration file to create a simulation instance
                sim_world = Simulation.from_file("config/config.yml")

                # Update the interaction coefficients based on the current iteration
                sim_world.container.neuron_factory.clocks_factory.proliferation_rate = prol_rate
                sim_world.container.neuron_factory.clocks_factory.differentiation_rate = diff_rate
                sim_world.container.neuron_factory.max_number_of_neurites = max_neurites

                # Run the simulation
                simulation(sim_world)

                # Plot and save the results (mark interactive as False to automatically  close the window)
                sim_world.container.animator.show(interactive=False)
                sim_world.container.animator.save_screenshot(f"output/{folder_name}/replicate{i}.png")
                sim_world.container.animator.plotter.close()


if __name__ == "__main__":
    main()