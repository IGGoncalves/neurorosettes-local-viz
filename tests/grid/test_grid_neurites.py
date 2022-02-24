"""Script to simulate the formation of a neuroblastoma rosette by cells with radial neurites (2D)"""
from typing import Dict

from neurorosettes.simulation import Container

# Time constants
STEP: float = 0.1
TOTAL_TIME: float = 20.0

# Domain constants
GRID_MIN: float = -60.0
GRID_MAX: float = 60.0
GRID_STEP: float = 20.0


def create_tissue(container: Container) -> None:
    # Populate environment with cells
    container.create_new_neuron(coordinates=[2.0, 10.0, 0.0],
                                outgrowth_axis=[1.0, 0.0, 0.0],
                                differentiation_grade=2)
    container.create_new_neuron(coordinates=[45.0, -12.0, 0.0])

    container.animator.set_camera(height=300.0)
    container.animator.show()


if __name__ == "__main__":
    # Initialize simulation objects
    sim_world = Container(grid_range=[GRID_MIN, GRID_MAX, GRID_STEP])
    # Create tissue
    create_tissue(sim_world)
    # Print neighbors
    print(sim_world.grid.get_close_objects(sim_world.neurons[0].neurites[-1].distal_point))
    # Plot the results (mark interactive as False to automatically  close the window)
    sim_world.animator.show(interactive=True)
