# Creating a model

## Basic model structure
Models can easily be created following the code structure presented below:

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
time:
  total_time: 100.0
  step: 0.1

domain:
  min: -60.0
  max: 60.0
  step: 20.0

use_2d: True

objects:
  cell_radius: 8.0
  cell_interaction_factor: 1.25
  neurite_radius: 0.5
  neurite_interaction_factor: 1.25
  neurite_spring_constant: 5.0
  neurite_default_length: 10.0

interactions:
  type: "potentials"
  sphere_sphere_adhesion: 4.0
  sphere_sphere_repulsion: 50.0
  sphere_sphere_smoothness: 1
  sphere_cylinder_adhesion: 5.0
  sphere_cylinder_repulsion: 10.0
  sphere_cylinder_smoothness: 1
  cylinder_cylinder_adhesion: 10.0
  cylinder_cylinder_repulsion: 10.0
  cylinder_cyliner_smoothness: 1
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
