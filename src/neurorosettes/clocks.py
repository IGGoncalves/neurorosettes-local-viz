"""This module deals with the proliferation, differentiation and death cycles"""
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


class Clock:
    """
    Class to represent a cell's internal clock.

    Cells start in an inactive status and can be flagged for an event based
    on a rate. Alternatively, this status can be manually induced.
    Active states can also be blocked. When a rate is not defined, cells remain
    in an inactive state.

    Parameters
    ----------
    rate
        The rate of the clock's event
    """

    def __init__(self, rate: Optional[float] = None) -> None:
        self.rate = rate
        self.signal = False
        self.block = False if rate else True

    def set_rate(self, rate: float) -> None:
        """
        Sets the clock's event rate.
        
        Parameters
        ----------
        rate
            The clock's event rate
        """
        self.rate = rate

    def advance_clock(self, timestep: float) -> None:
        """Updates the cell cycle status based on the proliferation rate (may happen or not)"""
        if self.block:
            return
        probability = timestep * self.rate
        self.signal = np.random.rand() <= probability

    def flag(self) -> None:
        """Updates the cell cycle to proliferation (will always happen)"""
        self.signal = True

    def remove_flag(self) -> None:
        """Updates the cell cycle to arrest (will always happen)"""
        self.signal = False

    def block(self) -> None:
        """Activates the cycle block to avoid proliferation"""
        self.signal = False
        self.block = True


class CellClocks:
    """
    Class to represent all of the internal biological clocks of a cell.

    Parameters
    ----------
    proliferation_rate
        A cell's proliferation rate.
    death_rate
        A cell's death rate.
    differentiation_rate
        A cell's differentiation rate.
    """

    def __init__(self, proliferation_rate, death_rate, differentiation_rate) -> None:
        self.cycle_clock = Clock(proliferation_rate)
        self.death_clock = Clock(death_rate)
        self.differentiation_clock = Clock(differentiation_rate)

    def set_proliferation_clock(self, proliferation_rate: float):
        """
        Sets the proliferation rate.

        Parameters
        ----------
        proliferation_rate
            A cell's proliferation rate
        """
        self.cycle_clock.set_rate(proliferation_rate)

    def set_death_clock(self, death_rate: float):
        """
        Sets the death rate.

        Parameters
        ----------
        death_rate
            A cell's proliferation rate
        """
        self.death_clock.set_rate(death_rate)

    def set_differentiation_clock(self, differentiation_rate: float):
        """
        Sets the differentiation rate.

        Parameters
        ----------
        differentiation_rate
            A cell's differentiation rate
        """
        self.differentiation_clock.set_rate(differentiation_rate)

    def set_clocks(self, proliferation_rate: float, death_rate: float, differentiation_rate: float) -> None:
        """
        Sets the rates for all biological events (proliferation, death and differentiation).

        Parameters
        ----------
        proliferation_rate
            A cell's proliferation rate
        death_rate
            A cell's proliferation rate
        differentiation_rate
            A cell's differentiation rate
        """
        self.set_proliferation_clock(proliferation_rate)
        self.set_differentiation_clock(differentiation_rate)
        self.set_death_clock(death_rate)

    def get_clock_rates(self) -> Tuple[float, float, float]:
        """Returns all the clock rates for biological events."""
        proliferation = self.cycle_clock.rate
        death = self.death_clock.rate
        differentiation = self.differentiation_clock.rate

        return proliferation, death, differentiation

    def advance_clocks(self, timestep: float):
        """
        Advance all biological clocks and flag cells for biological events.

        Parameters
        ----------
        timestep
            The simulation time step (time passed between status updates)
        """
        self.cycle_clock.advance_clock(timestep)
        self.death_clock.advance_clock(timestep)
        self.differentiation_clock.advance_clock(timestep)

    def block_all_clocks(self):
        """Blocks the cycling, death and differentiation clocks."""
        self.cycle_clock.block()
        self.death_clock.block()
        self.differentiation_clock.block()


@dataclass
class ClocksFactory:
    """
    Helper class to create instances of cell internal clocks based on the same rates.
    """
    proliferation_rate: float
    """The rate of proliferation of the cells."""
    death_rate: float
    """The rate of death of the cells."""
    differentiation_rate: float
    """The rate of differntiation of the cells."""

    def get_clocks(self) -> CellClocks:
        """Returns an instance of a CellClocks object with the chosen biology rates."""
        return CellClocks(self.proliferation_rate, 
                          self.death_rate, 
                          self.differentiation_rate)
