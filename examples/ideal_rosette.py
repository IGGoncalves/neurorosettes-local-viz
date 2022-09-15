"""Script to simulate the formation of Homer-Wright rosettes."""

from neurorosettes.simulation import Simulation, Container
import numpy as np


def create_tissue(container: Container) -> None:
    """Creates and registers new neurons in the simulation world."""
    radius = 24.0
    number_of_cells = 9
    t = np.linspace(0, 2*np.pi, number_of_cells, endpoint=False)

    x = radius * np.cos(t)
    y = radius * np.sin(t)

    for x_coord, y_coord in zip(x, y):
        neuron = container.create_new_neuron(coordinates=[x_coord, y_coord, 0])
        neuron.set_outgrowth_axis(np.subtract(np.zeros(3), neuron.cell.position))

    #neuron.clocks.block_all_clocks()
    #neuron.cell.sphere.c("red")

def main():
    # Initialize simulation objects
    sim_world = Simulation.from_file("config.yml")
    # Create initial configuration
    create_tissue(sim_world.container)
    # Plot the current state of the simulation
    sim_world.container.animator.set_camera(height=150.0)
    sim_world.container.animator.show()
    #sim_world.container.animator.save_screenshot("../output/ideal_case_start_small")
    #sim_world.save_meshes("ideal_case_fast_start_2")
    # Run the simulation to check if springs work

    sim_time = sim_world.timer.get_progress_bar()

    for t in sim_time.range():
        sim_world.container.advance_cycles(sim_world.timer.step)
        sim_world.container.kill()
        sim_world.container.differentiate()
        sim_world.container.divide()
        # Solve interactions and draw the new object positions
        sim_world.container.solve_mechanics(sim_world.timer.step)
        sim_world.container.update_drawings()

        # Update the simulation time on the simulation window
        if t % 10 == 0:
            sim_world.container.animator.update_clock(t)

        if t % 150 == 0:
            sim_world.save_meshes(f"../output/anim_ideal_t{t}")

        # Print time to the console as a progressbar
        sim_world.timer.current_time += sim_world.timer.step
        sim_time.print()

    # Plot the results (mark interactive as False to automatically  close the window)
    sim_world.container.animator.show(interactive=True)
    # Plot the results (mark interactive as False to automatically  close the window)
    #sim_world.container.animator.show(interactive=False)
    #sim_world.container.animator.save_screenshot("../output/ideal_case_small_2")
    sim_world.save_meshes("../output/anim_ideal_final")



if __name__ == "__main__":
    main()