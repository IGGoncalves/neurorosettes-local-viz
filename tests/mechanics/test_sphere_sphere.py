"""Script to test sphere-sphere interactions."""
import time
import click

from neurorosettes.simulation import Simulation, Container

def create_tissue(container: Container) -> None:
    """Creates and registers new neurons in the simulation world."""
    container.create_new_neuron(coordinates=[2.0, 0.0, 0.0])
    container.create_new_neuron(coordinates=[-2.0, 0.0, 0.0])

    # Plot the current state of the simulation
    container.animator.set_camera(height=70.0)
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
