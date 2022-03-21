"""Script to test the spring components of neurites."""
from neurorosettes.simulation import Simulation, Container
from neurorosettes.utilities import HexagonalTissue, RectangularTissue
from neurorosettes.grid import OneLevelDensityCheck

# Cell-cell interactions
cell_cell_adhesion = [5.0, 50.0, 100.0]
cell_cell_repulsion = [15.0, 50.0, 100.0]

# Initial cell positions
TISSUE = HexagonalTissue(size=160, spacing=16).get_coordinates()
# Contact inhibition function
DENSITY_CHECK = OneLevelDensityCheck(max_neighbors=19)

def set_density_check(container: Container) -> None:
    """Adds the contact inhibition function to the simulation."""
    container.set_density_check(DENSITY_CHECK)

def create_tissue(container: Container) -> None:
    """Creates and registers new neurons in the simulation world."""
    # Create a neuron with two neurites
    for cell in TISSUE:
        container.create_new_neuron(coordinates=cell)

def simulation(sim_world: Simulation) -> None:
    # Add additional functions to the simulation (i.g., contact inhibition)
    set_density_check(sim_world.container)
    # Create initial configuration
    create_tissue(sim_world.container)
    # Plot the current state of the simulation
    sim_world.container.animator.set_camera(height=600.0)
    sim_world.container.animator.show()
    # Run the simulation to check if springs work
    sim_world.run()
    # Plot the results (mark interactive as False to automatically  close the window)
    sim_world.container.animator.show(interactive=False)
    sim_world.container.animator.plotter.screenshot("THIS A TEST")
    sim_world.container.animator.plotter.close()


if __name__ == "__main__":
    sim_world = Simulation.from_file("config.yml")
    sim_world.container.sphere_int.adhesion_coefficient = cell_cell_adhesion[0]
    sim_world.container.sphere_int.repulsion_coefficient = cell_cell_repulsion[0]
    print(sim_world.container.sphere_int.adhesion_coefficient)
    simulation(sim_world)
