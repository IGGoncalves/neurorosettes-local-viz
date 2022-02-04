"""This module deals with the proliferation, differentiation and death cycles"""
from typing import Optional, Tuple

import numpy as np


class Cycle:
    """Represents the cell cycle (arrest/proliferation)"""
    proliferation_rate: Optional[float]
    division_signal: bool
    cycle_block: bool

    def __init__(self) -> None:
        self.proliferation_rate = None
        self.division_signal = False
        self.cycle_block = True

    def set_proliferation_rate(self, proliferation_rate: float):
        self.proliferation_rate = proliferation_rate
        self.cycle_block = False

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
    """Represents the cell death status (alive/dead)"""
    death_rate: Optional[float]
    death_signal: bool
    death_block: bool

    def __init__(self) -> None:
        self.death_rate = None
        self.death_signal = False
        self.death_block = True

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
        self.death_block = True


class Differentiation:
    differentiation_rate: Optional[float]
    differentiation_signal: bool
    differentiation_grade: int
    differentiation_block: bool

    def __init__(self) -> None:
        self.differentiation_rate = None
        self.differentiation_signal = False
        self.differentiation_grade = 0
        self.differentiation_block = True

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
        self.differentiation_block = True


class CellClocks:
    cycle_clock: Cycle
    death_clock: Death
    differentiation_clock: Differentiation

    def __init__(self) -> None:
        self.cycle_clock = Cycle()
        self.death_clock = Death()
        self.differentiation_clock = Differentiation()

    def set_proliferation_clock(self, proliferation_rate: float):
        self.cycle_clock.set_proliferation_rate(proliferation_rate)

    def set_death_clock(self, death_rate: float):
        self.death_clock.set_death_rate(death_rate)

    def set_differentiation_clock(self, differentiation_rate: float):
        self.differentiation_clock.set_differentiation_rate(differentiation_rate)

    def set_clocks(self, proliferation_rate: float, death_rate: float, differentiation_rate: float) -> None:
        self.set_proliferation_clock(proliferation_rate)
        self.set_differentiation_clock(differentiation_rate)
        self.set_death_clock(death_rate)

    def get_clock_rates(self) -> Tuple[float, float, float]:
        proliferation = self.cycle_clock.proliferation_rate
        death = self.death_clock.death_rate
        differentiation = self.differentiation_clock.differentiation_rate

        return proliferation, death, differentiation

    def advance_clocks(self, timestep: float):
        self.cycle_clock.advance_cycle_clock(timestep)
        self.death_clock.advance_death_clock(timestep)
        self.differentiation_clock.advance_differentiation_clock(timestep)

    def block_all_clocks(self):
        self.cycle_clock.block_proliferation()
        self.death_clock.block_death()
        self.differentiation_clock.block_differentiation()
