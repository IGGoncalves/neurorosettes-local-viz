# Contact inhibition functions

## One-level contact inhibition

The one-level contact inhibition function evaluates the number of neighbors in the current cell of the domain grid, as well as the cells surrounding it. This function evaluates if the number of cell bodies registered in these cells is larger than a user-defined value.

## Two-level contact inhibition

The two-level contact inhibition function consists of a first level that does the same as the one-level contact inhibition function, plus an additional evaluation. In this second level, the number of neighbors inside a user-defined radius of interest is evaluated, to assess if cells are fully surrounded by other cells or if there is some direction in which growth is possible.