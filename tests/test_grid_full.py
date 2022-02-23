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


def create_tissue(container: Container) -> Dict[int, str]:
    # Populate environment with cells
    neuron1 = container.create_new_neuron(coordinates=[2.0, 10.0, 0.0], color="green")
    neuron2 = container.create_new_neuron(coordinates=[25.0, -12.0, 0.0], color="orange")
    neuron3 = container.create_new_neuron(coordinates=[25.0, 32.0, 0.0], color="blue")
    neuron4 = container.create_new_neuron(coordinates=[-35.0, -22.0, 0.0], color="red")

    neuron_dict = {id(neuron.cell): color
                   for neuron, color in zip([neuron1, neuron2, neuron3, neuron4],
                                            ["green", "orange", "blue", "red"])}

    container.animator.set_camera(height=200.0)
    container.animator.show()

    return neuron_dict


def print_neighbors(neuron_index: int, neuron_dict: Dict[int, str]) -> None:
    print(f"The {neuron_dict[id(sim_world.neurons[neuron_index].cell)]} cell is surrounded by (including itself):")
    for neighbor in sim_world.grid.get_close_objects(sim_world.neurons[i].cell.position):
        print(f"The {neuron_dict[id(neighbor)]} cell")


if __name__ == "__main__":
    # Initialize simulation objects
    sim_world = Container(grid_range=[GRID_MIN, GRID_MAX, GRID_STEP])
    # Create initial configuration
    id_dict = create_tissue(sim_world)
    for i, neuron in enumerate(sim_world.neurons):
        print_neighbors(i, id_dict)
        print("-------------------")
    # Plot the results (mark interactive as False to automatically  close the window)
    sim_world.animator.show(interactive=True)
