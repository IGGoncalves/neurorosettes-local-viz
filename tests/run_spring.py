"""Script to test the spring components of neurites."""
import time

from neurorosettes.simulation import Simulation, Container

def create_tissue(container: Container) -> None:
    """Creates and registers new neurons in the simulation world."""
    neuron = container.create_new_neuron(coordinates=[-0.0, 10.0, 0.0],
                                         outgrowth_axis=[1.0, 0.0, 0.0])
                                         
    neuron.create_first_neurite(container.object_factory)
    neurite = neuron.neurites[0]
    neurite.create_neurite_representation(container.animator)

    container.grid.register_neurite(neurite)

    # Plot the current state of the simulation
    container.animator.set_camera(height=150.0)
    container.animator.show()
    time.sleep(1)


def update_tissue(container: Container) -> None:
    """Updates object positions inside the simulation world."""
    # Move the cell to test if the springs pull the parts together
    container.move_cell(container.neurons[0], [-10.0, 10.0, 0.0])
    # Plot the current state of the simulation
    container.update_drawings()
    time.sleep(1)


def main():
    # Initialize simulation objects
    sim_world = Simulation.from_file("config.yml")
    sim_world.timer.total_time = 25.0
    # Create initial configuration
    create_tissue(sim_world.container)
    # Move things around to test the springs
    update_tissue(sim_world.container)
    # Run the simulation to check if springs work
    sim_world.run()
    # Plot the results (mark interactive as False to automatically  close the window)
    sim_world.container.animator.show(interactive=True)


if __name__ == "__main__":
    main()
