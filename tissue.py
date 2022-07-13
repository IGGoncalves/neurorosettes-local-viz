"""Script to perform a parameter study on the cell mechanics parameters that regulate the tissue architecture."""
from typing import Dict

from neurorosettes.simulation import Simulation, Container
from neurorosettes.utilities import HexagonalTissue
from neurorosettes.grid import OneLevelDensityCheck


# Contact inhibition function
DENSITY_CHECK = OneLevelDensityCheck(max_neighbors=19)


def set_density_check(container: Container) -> None:
    """Adds the contact inhibition function to the simulation."""
    container.set_density_check(DENSITY_CHECK)


def create_tissue(container: Container) -> None:
    """Creates and registers new neurons in the simulation world."""
    # Initial cell positions
    tissue = HexagonalTissue(size=220, spacing=15).get_coordinates()
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
    sim_world.container.animator.save_screenshot(f"output/hex/low_start")
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
    # Load the data from the configuration file to create a simulation instance
    sim_world = Simulation.from_file("config/config.yml")

    # Run the simulation
    simulation(sim_world)

    # Plot and save the results (mark interactive as False to automatically  close the window)
    sim_world.container.animator.show(interactive=False)
    sim_world.container.animator.save_screenshot(f"output/hex/low_end")
    sim_world.container.animator.plotter.close()


if __name__ == "__main__":
    main()