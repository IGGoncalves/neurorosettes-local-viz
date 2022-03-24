"""Script to test the spring components of neurites."""
from typing import Dict
from pathlib import Path
import json
from itertools import count

from neurorosettes.simulation import Simulation, Container
from neurorosettes.utilities import HexagonalTissue, RectangularTissue
from neurorosettes.grid import OneLevelDensityCheck

# Cell-cell interactions
cell_cell_adhesion = 5.0
cell_cell_repulsion = [15.0, 25.0, 50.0]
# Cell-neurite interactions
cell_neurite_adhesion = 5.0
cell_neurite_repulsion = [50.0]
# Neurite-neurite interactions
neurite_neurite_adhesion = 5.0
neurite_neurite_repulsion = [50.0]

# Replicates for each condition
NUMBER_OF_REPLICATES = 3

# Contact inhibition function
DENSITY_CHECK = OneLevelDensityCheck(max_neighbors=19)


def set_density_check(container: Container) -> None:
    """Adds the contact inhibition function to the simulation."""
    container.set_density_check(DENSITY_CHECK)


def create_tissue(container: Container) -> None:
    """Creates and registers new neurons in the simulation world."""
    # Initial cell positions
    tissue = HexagonalTissue(size=160, spacing=20).get_coordinates()
    # Create a neuron with two neurites
    for cell in tissue:
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


def create_params_dict(
    cca: float, ccr: float, cna: float, cnr: float, nna: float, nnr: float
) -> Dict[str, float]:
    return {
        "cell-cell adhesion": cca,
        "cell-cell repulsion": ccr,
        "cell-neurite adhesion": cna,
        "cell_neurite repulsion": cnr,
        "neurite-neurite adhesion": nna,
        "neurite-neurite repulsion": nnr,
    }


if __name__ == "__main__":
    sim_counter = count(start=24)
    for nnr in neurite_neurite_repulsion:
        for cnr in cell_neurite_repulsion:
            for ccr in cell_cell_repulsion:
                parameters = create_params_dict(
                    cell_cell_adhesion,
                    ccr,
                    cell_neurite_adhesion,
                    cnr,
                    neurite_neurite_adhesion,
                    nnr,
                )

                folder_name = next(sim_counter)
                filepath = Path(f"output/{folder_name}/params.json")
                filepath.parent.mkdir(parents=True, exist_ok=True)

                with open(filepath, "w") as outfile:
                    json.dump(parameters, outfile)

                for i in range(NUMBER_OF_REPLICATES):
                    sim_world = Simulation.from_file("config.yml")
                    sim_world.container.sphere_int.adhesion_coefficient = cell_cell_adhesion
                    sim_world.container.sphere_int.repulsion_coefficient = ccr
                    sim_world.container.sphere_cylinder_int.adhesion_coefficient = cell_neurite_adhesion
                    sim_world.container.sphere_cylinder_int.repulsion_coefficient = cnr
                    sim_world.container.cylinder_int.adhesion_coefficient = neurite_neurite_adhesion
                    sim_world.container.cylinder_int.repulsion_coefficient = nnr
                    simulation(sim_world)
                    # Plot the results (mark interactive as False to automatically  close the window)
                    sim_world.container.animator.show(interactive=False)
                    sim_world.container.animator.save_screenshot(f"output/{folder_name}/replicate{i}.png")
                    sim_world.container.animator.plotter.close()