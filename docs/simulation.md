# Creating a model

Models are defined by two main components: the `config.yaml` file and the main Python script. The configuration file contains all the basic model parameters required to run a standard model. The main Python script creates a simulation object based on the configuration file, deals with populating the simulation domain with neurons, and runs the simulation.

## Basic model structure

A standard script for a `neurorosettes` model can be divided in three groups:

- **Simulation initialization**: the imualtion parameters are defined and a `Simulation` object is created;

- **Tissue creation**: neurons are created and placed in the simulation domain;

- **Running the simulation**: the model is run for the time defined by the user. At every time point, the simulation objects' status are updated to model proliferation, differentiation and death, and their position is also updated based on the mechanical interactions between objects.

Each of these steps is described in further detail in the following sections.

Takking this into account, models can easily be created following the code structure presented below:

```python
"""Script to test the spring components of neurites."""
from neurorosettes.simulation import Simulation, Container


def create_tissue(container: Container) -> None:
    """Creates and registers new neurons in the simulation world."""
    pass

def main() -> None:
     # Initialize simulation objects
    sim_world = Simulation.from_file("config.yml")
    # Create initial configuration
    create_tissue(sim_world.container)
    # Run the simulation to check if springs work
    sim_world.run()
    # Plot the results (mark interactive as False to automatically  close the window)
    sim_world.container.animator.show(interactive=True)

if __name__ == "__main__":
    main()
```

### Loading data from the configuration file

All the data required to create a new simulation is stored in YAML configuration files.
This makes it easier for users to change parameters in a simulation without having to
interact with actual code, as the `Simulation` class can be initialized through the 
`from_file()` method, which takes the path to the configuration file as input.

The structure of a simulation config file is the following:

```yaml
# Simulation time [minutes]
time:
  total_time: 1000.0
  step: 0.1

# Domain options
domain:
  # Domain size [microns]
  boundaries:
    min: -300.0
    max: 300.0
    step: 20.0

  # 2D status [True for 2D / False for 3D]
  use_2d: True

  # Drag coefficent [Pa*s]
  drag_coefficient: 10.0

# Neuron objects
neurons:
  # Physical objects properties
  objects:
    cell_radius: 8.0
    cell_interaction_factor: 1.25
    neurite_radius: 0.5
    neurite_interaction_factor: 5.0
    neurite_spring_constant: 5.0
    neurite_default_length: 10.0

  # Cycling and differentitation rates [1/min]
  clocks:
    proliferation_rate: 0.0007
    death_rate: 0.00002
    differentiation_rate: 0.0006

  # Number of neurites allowed until differentiation is blocked 
  max_number_of_neurites: 3

# Object interactions
interactions:
  # Interactions may be based on simple adhesion/repulsion functions ["simple"]
  # or potentials functions ["potentials"]
  type: "potentials"
  sphere_sphere_adhesion: 5.0
  sphere_sphere_repulsion: 30.0
  # Smoothness factor should be commented for simple adhesion/repulsion
  sphere_sphere_smoothness: 1
  sphere_cylinder_adhesion: 5.0
  sphere_cylinder_repulsion: 30.0
  # Smoothness factor should be commented for simple adhesion/repulsion
  sphere_cylinder_smoothness: 1
  cylinder_cylinder_adhesion: 10.0
  cylinder_cylinder_repulsion: 30.0
    # Smoothness factor should be commented for simple adhesion/repulsion
  cylinder_cylinder_smoothness: 1
```

### Creating the tissue structure

Users should define how tissues are initialized by creating and placing neurons in the
simulation domain. This step of the process is not fixed, as users should have the
flexibility to create new tissues based on their study needs. However, some 
[initialization functions](tissue_initialization.md) are provided with the
package, to create structured grid tissues.

### Running the simulation

Simulations can be run by calling the `.run()` method for the `Simulation` instance.
This function will run the simulation for the total time specified in the `config.yml`
file. At each time point, the biological clocks are advanced, and then neurons are checked
for cell, death and proliferation. Mechanical interactions are also solved and object
positions are updated accordingly.

## Model customization

It is possible to add more functionalities to the model. For instance, users can choose to include [contact inhibition functions](contact_inhibition.md) that will check the number of neighbors of a cell before proliferation.
