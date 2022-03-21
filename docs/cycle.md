# Cell cycles

Neurons can have associated biological clocks that define when a neuron may proliferate,
differentiate or day. Each of these processes has a defined rate, and the probability
of a cell to undergo that pathway at a given point in time is defined by the equation
below:

$ P_{\textrm{process}} \approx r_{\textrm{process}} \Delta t$

At each time point, neurons clocks are advanced and neurons may be flagged for one of these
processes, according to the previously defined probabilities. Once a neuron has been flagged
for a process, a new object may be created, if we are dealing with proliferation or
differentiation, or the current neuron may be removed, if it has been flagged for death.

Cells may have their biological events blocked by the user or in response to external conditions, by calling
the `block()` method. Users may also flag a cell for a given event, overriding
the proliferation rate, using the `flag()` method. [Contact inhibition
functions](contact_inhibition.md) may be used to block proliferation when the neuron density
is too high to enable division.

## Proliferation

Proliferation is possible for neurons with a cell body and no neurites. When a neuron
divides, a daughter neuron is initialized and registered in the simulation domain. This
neuron inherits the mother's mechanical properties and cycle rates.

## Death

Neuron death can occur regardless of the number of neurites that a neuron may have. When
a neuron is flagged for death, the cell body and neurites are removed from the simulation
domain. 

## Differentiation

Differentiation is characterized by the formation of a primary neurite or the creation
of a secondary neurite, that extends on a previous one. When it occurs, proliferation is blocked
for the neuron, as it has been commited for differentiation.
