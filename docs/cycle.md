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

## Proliferation

Proliferation is possible for neurons with a cell body and no neurites. When a neuron
divides, a daughter neuron is initialized and registered in the simulation domain. This
neuron inherits the mother's mechanical properties and cycle rates. Cells may have their
proliferation cycle blocked by the user or in response to external conditions, by calling
the `block_proliferation()` method. Users may also flag a cell for proliferation, overriding
the proliferation rate, using the `flag_for_proliferation()` method. [Contact inhibition
functions](contact_inhibition.md) may be used to block proliferation when the neuron density
is too high to enable division.

## Death

Neuron death can occur regardless of the number of neurites that a neuron may have. When
a neuron is flagged for death, the cell body and neurites are removed from the simulation
domain. As for proliferation, the `block_death()` and `flag_for_death` methods are available.

## Differentiation

Differentiation is characterized by the formation of a primary neurite or the creation
of a secondary neurite, that extends on a previous one. As before, the `block_differentiation()` 
and `flag_for_differentiation` methods are available.

To avoid the overlap between a new
neurite and an existing one, which is more difficult to solve than the overlap between
too cell bodies, new neurites are only formed if they do not overlap. If there is overlap,
the neurite is shortened to avoid the intersection between the new neurites. If this is not
possible either, a new axis of growth is randomly selected.