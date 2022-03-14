"""Script to test the spring components of neurites."""
from neurorosettes.simulation import Simulation, Container
from neurorosettes.utilities import RectangularTissue
from neurorosettes.grid import OneLevelDensityCheck


TISSUE = RectangularTissue(size=160, spacing=20).get_coordinates()
DENSITY_CHECK = OneLevelDensityCheck(max_neighbors=19)

def set_density_check(container: Container) -> None:
    container.set_density_check(DENSITY_CHECK)

def create_tissue(container: Container) -> None:
    """Creates and registers new neurons in the simulation world."""
    # Create a neuron with two neurites
    for cell in TISSUE:
        container.create_new_neuron(coordinates=cell)

    # Plot the current state of the simulation
    container.animator.set_camera(height=600.0)
    container.animator.show()


if __name__ == "__main__":
    # Initialize simulation objects
    sim_world = Simulation.from_file("config.yml")
    set_density_check(sim_world.container)
    # Create initial configuration
    create_tissue(sim_world.container)
    # Run the simulation to check if springs work
    sim_world.run()
    # Plot the results (mark interactive as False to automatically  close the window)
    sim_world.container.animator.show(interactive=True)
