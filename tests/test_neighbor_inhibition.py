"""Script to simulate the formation of a neuroblastoma rosette by cells with radial neurites (2D)"""
from typing import Tuple
import time

import numpy as np
import vedo

from neurorosettes.neurons import Neuron
from neurorosettes.subcellular import CellBody
from neurorosettes.simulation import Container
from neurorosettes import utilities


def create_tissue() -> list:
    """Returns a list of randomly distributed positions for a given number of cells"""
    coordinates = [np.array([x + 8 * (i % 2), y, 0])
                   for x in range(-25, 5, 25)
                   for i, y in enumerate(range(-50, 51, 25))]
    return coordinates


def create_legend_item(status: str, color: str, pos: Tuple[int]) -> Tuple[vedo.Circle, vedo.Text2D]:
    label = vedo.Circle(pos=pos, r=5, c=color, alpha=1, res=120)
    text = vedo.Text2D(status, pos=(0.83, 0.66), c="black")

    return label, text


# Define time variables
timestep = 0.1
total_time = 300
pb = utilities.get_simulation_timer(total_time, timestep)

clock = vedo.Text2D("Simulation step: 0", pos="top right", c="black", font="Courier")
legend = vedo.shapes.Rectangle(p1=(70, 50, 0), p2=(122, -30), c='gray7', alpha=1)
legend_outline = vedo.shapes.Rectangle(p1=(69.5, 50.5, 0), p2=(122.5, -30.5, 0), c='black', alpha=1)

cycling_legend = vedo.shapes.Circle(pos=(80, 40, 0), r=5, c='red', alpha=1, res=120)
cycling_text = vedo.Text2D("Live", pos=(0.83, 0.66), c="black", font="Courier")

blocked_legend = vedo.shapes.Circle(pos=(80, 25, 0), r=5, c='black', alpha=1, res=120)
blocked_text = vedo.Text2D("Blocked", pos=(0.83, 0.61), c="black", font="Courier")

active_legend = vedo.shapes.Circle(pos=(80, 10, 0), r=5, c='blue', alpha=1, res=120)
active_text = vedo.Text2D("Active", pos=(0.83, 0.55), c="black", font="Courier")

neighbor_legend = vedo.shapes.Circle(pos=(80, -5, 0), r=5, c='yellow', alpha=1, res=120)
neighbor_text = vedo.Text2D("Neighbor", pos=(0.83, 0.49), c="black", font="Courier")

nearby_legend = vedo.shapes.Circle(pos=(80, -20, 0), r=5, c='brown', alpha=1, res=120)
nearby_text = vedo.Text2D("Neighbor+", pos=(0.83, 0.44), c="black", font="Courier")

# Initialize simulation objects
container = Container(timestep=timestep,
                      simulation_2d=True,
                      viscosity=7.96,
                      grid=[-160, 160, 20])

for position in create_tissue():
    # Populate environment with cells
    neuron = Neuron()
    neuron.create_cell(coordinates=np.array(position))
    neuron.set_outgrowth_axis(utilities.get_random_unit_vector(two_dimensions=True))
    neuron.clocks.set_clocks(0.04, 0.00000001, 0.0000002)
    container.register_neuron(neuron)

container.animator.plotter += clock
container.animator.plotter += legend_outline
container.animator.plotter += legend
container.animator.plotter += cycling_legend
container.animator.plotter += cycling_text
container.animator.plotter += blocked_legend
container.animator.plotter += blocked_text
container.animator.plotter += active_legend
container.animator.plotter += active_text
container.animator.plotter += neighbor_legend
container.animator.plotter += neighbor_text
container.animator.plotter += nearby_legend
container.animator.plotter += nearby_text

container.animator.plotter.show(resetcam=False, interactive=False)
blocked_cells = []

# Run and plot simulation
for t in pb.range():
    if t % 5 == 0:
        container.advance_cycles()
        container.kill()
        container.differentiate()

        for neuron in container.neurons:
            if not neuron.ready_for_division:
                continue

            cell_neighbors = [neighbor
                              for neighbor in container.grid.get_close_objects(neuron.cell.position)
                              if isinstance(neighbor, CellBody)]

            if len(cell_neighbors) > 15:

                acceptable_radius = neuron.cell.mechanics.radius * 2.1
                nearby = [neighbor for neighbor in cell_neighbors
                          if np.linalg.norm(neighbor.position - neuron.cell.position) <= acceptable_radius]

                if len(nearby) < 7:
                    continue

                neuron.cell.sphere.c("blue")
                for neighbor in cell_neighbors:
                    if neighbor is neuron.cell:
                        continue
                    if neighbor in nearby:
                        neighbor.sphere.c("brown")
                    else:
                        neighbor.sphere.c("yellow")

                container.update_drawings()
                time.sleep(2)
                for neighbor in cell_neighbors:
                    if neighbor is neuron.cell:
                        neighbor.sphere.c("black")
                        blocked_cells.append(neighbor)
                    else:
                        if neighbor not in blocked_cells:
                            neighbor.sphere.c("red")
                        else:
                            neighbor.sphere.c("black")

                container.update_drawings()
                neuron.clocks.cycle_clock.division_signal = False
                neuron.clocks.cycle_clock.cycle_block = True

        container.divide()
    container.update_cell_positions()
    container.update_drawings()

    if t % 50 == 0:
        clock.text(f"Simulation step: {t}")

    pb.print()

container.animator.plotter.show(interactive=True)
