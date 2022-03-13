"""Script to test the spring components of neurites."""
import time

from neurorosettes.simulation import Simulation, Container
from neurorosettes.utilities import RectangularTissue


tissue = RectangularTissue(160, spacing=20).get_coordinates()

def create_tissue(container: Container) -> None:
    """Creates and registers new neurons in the simulation world."""
    # Create a neuron with two neurites
    for cell in tissue:
        container.create_new_neuron(coordinates=cell)


    # Plot the current state of the simulation
    container.animator.set_camera(height=600.0)
    container.animator.show()
    time.sleep(1)



if __name__ == "__main__":
    # Initialize simulation objects
    sim_world = Simulation.from_file("config.yml")
    # Create initial configuration
    create_tissue(sim_world.container)
    # Run the simulation to check if springs work
    sim_world.run()
    # Plot the results (mark interactive as False to automatically  close the window)
    sim_world.container.animator.show(interactive=True)
