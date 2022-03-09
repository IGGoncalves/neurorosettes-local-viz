import unittest

import numpy as np

from neurorosettes.grid import UniformGrid
from neurorosettes.neurons import ObjectFactory


factory = ObjectFactory(cell_radius=8.0, cell_interaction_factor=1.25, neurite_radius=0.5,
                        neurite_interaction_factor=1.25, neurite_spring_constant=10.0, neurite_default_length=10.0)

class GridTest(unittest.TestCase):
    def test_interpolate_idx(self):
        # Create a 2x2 grid with value bounds [-20, 0, (20)]
        grid = UniformGrid(min=-20.0, max=20.0, step=20.0)
        # Test a position on the top right cell
        i, j, _ = grid.interpolate_idx([5.0, 15.0, 0.0])
        self.assertEqual((i, j), (1, 1))
        # Test a position on the bottom left cell
        i, j, _ = grid.interpolate_idx([-5.0, -5.0, 0.0])
        self.assertEqual((i, j), (0, 0))

    def test_register_cell(self):
        grid = UniformGrid(min=-20.0, max=20.0, step=20.0)
        cell = factory.get_cell_body([-5.0, 5.0, 0.0])
        grid.register_cell(cell)

        self.assertEqual(grid.get_objects_in_voxel(0, 1), [cell])


if __name__ == '__main__':
    unittest.main()