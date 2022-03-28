# Running examples

We offer some scripts that can be run out of the box once the `neurorosettes` package has been installed.

- **rosette_formation** - Runs a single simulation to reproduce the formation of rosettes in tissues nervous system.
- **sweep_mechanics** - Runs a parameter study on the adhesion and repulsion coefficients that regulate physical interactions between cell bodies and neurites. Can be adapted to perform parameter studies on other variables.

Additional scripts to evaluate mechanical interactions:

- **cell_cell** - Runs a simulation with two cells that start at a close distance and are moved based on the interaction forces acting on them.
- **cell-neurite** - Similar to `cell_cell` but with a cell and a neurite.
- ** neurite-neurite** - Similar to `cell_cell` but with two neurites.

## Running examples from the commandline

All of the scripts mentioned in this seciton work in a similar way. To run them, their name should be called in the commandline. It is expected that they are run from a location that contains a folder name `config` which should
have a configuration YAML file inside. If this is not the case, the script accepts an additional argument `config_path` that should point to where this file is stored. For instance, if we wanted to run the cell_cell script and
point to a config file in the same directory, we would run:

```sh
cell_cell --config_path="config.yml"
```