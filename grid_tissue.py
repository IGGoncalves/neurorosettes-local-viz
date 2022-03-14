"""Script to test the spring components of neurites."""
from neurorosettes.simulation import Simulation, Container
from neurorosettes.utilities import HexagonalTissue
from neurorosettes.grid import OneLevelDensityCheck


TISSUE = HexagonalTissue(size=160, spacing=18).get_coordinates()
DENSITY_CHECK = OneLevelDensityCheck(max_neighbors=19)
DRAG_COEFFICIENT = 10.0

def set_density_check(container: Container) -> None:
    container.set_density_check(DENSITY_CHECK)

def set_drag_coefficient(container: Container) -> None:
    container.drag_coefficient = DRAG_COEFFICIENT

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
    set_drag_coefficient(sim_world.container)
    # Create initial configuration
    create_tissue(sim_world.container)
    # Run the simulation to check if springs work
    sim_world.run()
    # Plot the results (mark interactive as False to automatically  close the window)
    sim_world.container.animator.show(interactive=True)
