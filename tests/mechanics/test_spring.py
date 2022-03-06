"""Script to test the spring components of neurites."""
import time

from neurorosettes.simulation import Simulation, Container


def create_tissue(container: Container) -> None:
    """Creates and registers new neurons in the simulation world."""
    # Create a neuron with two neurites
    container.create_new_neuron(coordinates=[-0.0, 10.0, 0.0],
                                outgrowth_axis=[1.0, 0.0, 0.0],
                                differentiation_grade=2)

    # Create a neuron with one neurite
    container.create_new_neuron(coordinates=[-5.0, -10.0, 0.0],
                                outgrowth_axis=[1.0, 0.0, 0.0],
                                differentiation_grade=1)

    # Plot the current state of the simulation
    container.animator.set_camera(height=150.0)
    container.animator.show()
    time.sleep(1)


def update_tissue(container: Container) -> None:
    """Updates object positions inside the simulation world."""
    # Move the cell to test if the springs pull the parts together
    container.move_cell(container.neurons[0], [-10.0, 10.0, 0.0])
    # Move the cell to test if the spring pushes the parts away
    container.move_cell(container.neurons[1], [2.0, -10.0, 0.0])
    # Plot the current state of the simulation
    container.update_drawings()
    time.sleep(2)


def run_simulation(simulation: Simulation) -> None:
    """Runs the entire simulation by solving the mechanics at each time point."""
    sim_time = simulation.timer.get_progress_bar()

    for t in sim_time.range():
        # Solve interactions and draw the new object positions
        simulation.container.solve_mechanics(simulation.timer.step)
        simulation.container.update_drawings()

        # Update the simulation time on the simulation window
        if t % 10 == 0:
            simulation.container.animator.update_clock(t)

        # Print time to the console as a progressbar
        sim_time.print()


if __name__ == "__main__":
    # Initialize simulation objects
    sim_world = Simulation.from_file("config.yml")
    # Create initial configuration
    create_tissue(sim_world.container)
    # Move things around to test the springs
    update_tissue(sim_world.container)
    # Run the simulation to check if springs work
    run_simulation(sim_world)
    # Plot the results (mark interactive as False to automatically  close the window)
    sim_world.container.animator.show(interactive=True)
