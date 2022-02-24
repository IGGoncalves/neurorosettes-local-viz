"""Script to test the spring components of neurites."""
import time

from neurorosettes.simulation import Container
from neurorosettes.utilities import get_simulation_timer

# Time constants
STEP: float = 0.1
TOTAL_TIME: float = 20.0

# Domain constants
GRID_MIN: float = -60.0
GRID_MAX: float = 60.0
GRID_STEP: float = 20.0


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


def run_simulation(time_step: float, total_time: float, container: Container) -> None:
    """Runs the entire simulation by solving the mechanics at each time point."""
    sim_time = get_simulation_timer(total_time, time_step)

    for t in sim_time.range():
        # Solve interactions and draw the new object positions
        container.solve_mechanics(time_step)
        container.update_drawings()

        # Update the simulation time on the simulation window
        if t % 10 == 0:
            container.animator.update_clock(t)

        # Print time to the console as a progressbar
        sim_time.print()


if __name__ == "__main__":
    # Initialize simulation objects
    sim_world = Container(grid_range=[GRID_MIN, GRID_MAX, GRID_STEP])
    # Create initial configuration
    create_tissue(sim_world)
    # Move things around to test the springs
    update_tissue(sim_world)
    # Run the simulation to check if springs work
    run_simulation(STEP, TOTAL_TIME, sim_world)
    # Plot the results (mark interactive as False to automatically  close the window)
    sim_world.animator.show(interactive=True)
