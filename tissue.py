"""Script to perform a parameter study on the cell mechanics parameters that regulate the tissue architecture."""

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
    tissue = HexagonalTissue(size=190, spacing=20).get_coordinates()
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
    sim_world.container.animator.set_camera(height=400.0)
    sim_world.container.animator.show()
    #sim_world.container.animator.save_screenshot(f"output/start")
    #sim_world.save_meshes("output/start")
    # Run the simulation to check if springs work
    sim_world.run()

def main():
    # Load the data from the configuration file to create a simulation instance
    sim_world = Simulation.from_file("config/config.yml")

    # Run the simulation
    simulation(sim_world)

    # Plot and save the results (mark interactive as False to automatically close the window)
    #sim_world.container.animator.show(interactive=True)
    sim_world.container.animator.save_screenshot(f"output/no_diff")
    sim_world.save_meshes("output/no_diff")
    #sim_world.container.animator.plotter.close()


if __name__ == "__main__":
    main()