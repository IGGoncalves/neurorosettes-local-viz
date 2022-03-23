"""Script to test the spring components of neurites."""
from neurorosettes.simulation import Simulation, Container
from neurorosettes.utilities import HexagonalTissue, RectangularTissue
from neurorosettes.grid import OneLevelDensityCheck

# Cell-cell interactions
cell_cell_adhesion = 5.0
cell_cell_repulsion = [15.0, 25.0, 50.0]
# Cell-neurite interactions
cell_neurite_adhesion = 5.0
cell_neurite_repulsion = [15.0, 25.0, 50.0]
# Neurite-neurite interactions
neurite_neurite_adhesion = 5.0
neurite_neurite_repulsion = [15.0, 25.0, 50.0]

# Initial cell positions
TISSUE = HexagonalTissue(size=160, spacing=20).get_coordinates()
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


if __name__ == "__main__":
    for ccr in cell_cell_repulsion:
        for i in range(3):
            sim_world = Simulation.from_file("config.yml")
            sim_world.container.sphere_int.adhesion_coefficient = cell_cell_adhesion
            sim_world.container.sphere_int.repulsion_coefficient = ccr
            simulation(sim_world)
            # Plot the results (mark interactive as False to automatically  close the window)
            sim_world.container.animator.show(interactive=False)
            sim_world.container.animator.save_screenshot(f"{cell_cell_adhesion}_{ccr}_{i}.png")
            sim_world.container.animator.plotter.close()
    