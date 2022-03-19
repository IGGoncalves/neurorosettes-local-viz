"""This module deals with the proliferation, differentiation and death cycles"""
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


class Cycle:
    """
    Class to represent a cell's internal cycling clock.

    Cells start in a non-proliferative status and can be flagged for proliferation
    based on a proliferation rate. Alternatively, this status can be manually induced.
    Proliferation can also be blocked.

    Parameters
    ----------
    proliferation_rate
        The rate of proliferation of a cell
    """

    def __init__(self, proliferation_rate: Optional[float] = None) -> None:
        self.proliferation_rate = proliferation_rate
        self.division_signal = False
        self.cycle_block = False if proliferation_rate else True

    def advance_cycle_clock(self, timestep: float) -> None:
        """Updates the cell cycle status based on the proliferation rate (may happen or not)"""
        if self.cycle_block:
            return
        division_probability = timestep * self.proliferation_rate
        self.division_signal = np.random.rand() <= division_probability

    def flag_for_proliferation(self) -> None:
        """Updates the cell cycle to proliferation (will always happen)"""
        self.division_signal = True

    def remove_proliferation_flag(self) -> None:
        """Updates the cell cycle to arrest (will always happen)"""
        self.division_signal = False

    def block_proliferation(self) -> None:
        """Activates the cycle block to avoid proliferation"""
        self.division_signal = False
        self.cycle_block = True


class Death:
    """
    Class to represent a cell's internal death clock.

    Cells start in an alive status and can be flagged for death
    based on a death rate. Alternatively, this status can be manually induced.
    Death can also be blocked.

    Parameters
    ----------
    death_rate
        The rate of death of a cell
    """

    def __init__(self, death_rate: Optional[float] = None) -> None:
        self.death_rate = death_rate
        self.death_signal = False
        self.death_block = False if death_rate else True

    def set_death_rate(self, death_rate: float) -> None:
        self.death_rate = death_rate
        self.death_block = False

    def advance_death_clock(self, timestep: float) -> None:
        """Updates the cell status based on the death rate (may happen or not)"""
        if self.death_signal or self.death_block:
            return
        death_probability = timestep * self.death_rate
        self.death_signal = np.random.rand() <= death_probability

    def flag_for_death(self) -> None:
        """Updates the cell status to death (will always happen)"""
        self.death_signal = True

    def block_death(self) -> None:
        """Sets the death block status to True, to avoid cell death"""
        self.death_block = True


class Differentiation:
    """
    Class to represent a cell's internal differentiation clock.

    Cells start in an undifferentiated status and can be flagged for differentiation
    based on a differentiation rate. Alternatively, this status can be manually induced.
    Differentiation can also be blocked.

    Parameters
    ----------
    differentiation_rate
        The rate of death of a cell
    """

    def __init__(self, differentiation_rate: Optional[float] = None) -> None:
        self.differentiation_rate = differentiation_rate
        self.differentiation_signal = False
        self.differentiation_grade = 0
        self.differentiation_block = False if differentiation_rate else True

    def set_differentiation_rate(self, differentiation_rate: float):
        self.differentiation_rate = differentiation_rate
        self.differentiation_block = False

    def advance_differentiation_clock(self, timestep: float) -> None:
        """Updates the cell status based on the death rate (may happen or not)"""
        if self.differentiation_block:
            return
        death_probability = timestep * self.differentiation_rate
        self.differentiation_signal = np.random.rand() <= death_probability

    def increase_differentiation_grade(self):
        """Updates the differentiation grade and resets the signal to allow differentiation"""
        if not self.differentiation_signal:
            return
        self.differentiation_grade += 1
        self.differentiation_signal = False

    def block_differentiation(self) -> None:
        """Sets the differentiation block status to True, to avoid cell differntiation"""
        self.differentiation_block = True


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
        self.cycle_clock = Cycle(proliferation_rate)
        self.death_clock = Death(death_rate)
        self.differentiation_clock = Differentiation(differentiation_rate)

    def set_proliferation_clock(self, proliferation_rate: float):
        """
        Sets the proliferation rate.

        Parameters
        ----------
        proliferation_rate
            A cell's proliferation rate
        """
        self.cycle_clock.set_proliferation_rate(proliferation_rate)

    def set_death_clock(self, death_rate: float):
        """
        Sets the death rate.

        Parameters
        ----------
        death_rate
            A cell's proliferation rate
        """
        self.death_clock.set_death_rate(death_rate)

    def set_differentiation_clock(self, differentiation_rate: float):
        """
        Sets the differentiation rate.

        Parameters
        ----------
        differentiation_rate
            A cell's differentiation rate
        """
        self.differentiation_clock.set_differentiation_rate(differentiation_rate)

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
        proliferation = self.cycle_clock.proliferation_rate
        death = self.death_clock.death_rate
        differentiation = self.differentiation_clock.differentiation_rate

        return proliferation, death, differentiation

    def advance_clocks(self, timestep: float):
        """
        Advance all biological clocks and flag cells for biological events.

        Parameters
        ----------
        timestep
            The simulation time step (time passed between status updates)
        """
        self.cycle_clock.advance_cycle_clock(timestep)
        self.death_clock.advance_death_clock(timestep)
        self.differentiation_clock.advance_differentiation_clock(timestep)

    def block_all_clocks(self):
        """Blocks the cycling, death and differentiation clocks."""
        self.cycle_clock.block_proliferation()
        self.death_clock.block_death()
        self.differentiation_clock.block_differentiation()


@dataclass
class ClocksFactory:
    """
    Class to create instances of cell internal clocks based on the same rates.
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
