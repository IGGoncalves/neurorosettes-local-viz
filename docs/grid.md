# Domain discretization grid

To avoid computing interactions with a large number of neighbors, the simulation container
can be defined to have a `Grid` object, that discretizes the domain. Subcellular components
(cell bodies and neurites) are registered to the grid, which is updated at every time step
based on the new positions of the elements inside the grid. Objects can be retrieved from
the grid based on the index of a cell of the grid or the corresponding position.

## UniformGrid class

The current implementation of the `Grid` class is based on a uniform grid composed of cells
that have the same length on every dimension. The grid can be defined in 2D or 3D, and it must
be initialized with the bounds of the domain and the step size, which is the length of the
grid cells.

```python
from neurorosettes.grid import UniformGrid

grid = UniformGrid(min=-60.0, max=60.0, step=20.0)
```

