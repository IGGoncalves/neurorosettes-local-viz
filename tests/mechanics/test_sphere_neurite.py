"""Script to test sphere-neurite interactions"""
import time
import click

from neurorosettes.simulation import Simulation, Container


def create_tissue(container: Container) -> None:
    """Creates and registers new neurons in the simulation world."""
    neuron = container.create_new_neuron(coordinates=[0.0, -10.0, 0.0],
                                         outgrowth_axis=[0.0, 1.0, 0.0])

    neuron.create_first_neurite(container.object_factory)
    neuron.create_secondary_neurite(container.object_factory)
    for neurite in neuron.neurites:
        neurite.create_neurite_representation(container.animator)
        container.grid.register_neurite(neurite)

    container.create_new_neuron(coordinates=[4.0, 15.0, 0.0])

    # Plot the current state of the simulation
    container.animator.set_camera(height=100.0)
    container.animator.show()
    time.sleep(1)


@click.command()
@click.option("--config_path", default="config.yml", help="Configuration file path.")
def main(config_path):
    # Initialize simulation objects
    sim_world = Simulation.from_file(config_path)
    sim_world.timer.total_time = 25.0
    # Create initial configuration
    create_tissue(sim_world.container)
    # Run the simulation to check if springs work
    sim_world.run()
    # Plot the results (mark interactive as False to automatically  close the window)
    sim_world.container.animator.show(interactive=True)


if __name__ == "__main__":
    main()
