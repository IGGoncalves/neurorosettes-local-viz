import yaml
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, validator


class TimerValidator(BaseModel):
    """Class to parse and validate data on simulation time."""

    total_time: float
    step: float

    @validator("total_time", "step")
    def not_negative(cls, v):
        """Validates that simulation time data is not negative."""
        if v < 0.0:
            raise ValueError("Time data should be positive.")
        return v


class DomainValidator(BaseModel):
    """Class to parse and validate data on simulation space."""

    min: float
    max: float
    step: float

    @validator("step")
    def not_negative(cls, v):
        """Validates that simulation space data is not negative."""
        if v < 0.0:
            raise ValueError("Domain step should be positive.")
        return v


class ObjectValidator(BaseModel):
    """Class to parse and validate data on physical objects."""

    cell_radius: float
    cell_interaction_factor: float
    neurite_radius: float
    neurite_interaction_factor: float
    neurite_spring_constant: float
    neurite_default_length: float

    @validator("cell_radius")
    def not_negative(cls, v):
        """Validates that physical objects data is not negative."""
        if v < 0.0:
            raise ValueError("Physical object data should be positive.")
        return v


class ClocksValidator(BaseModel):
    """Class to parse and validate data on biological clocks."""

    proliferation_rate: float
    death_rate: float
    differentiation_rate: float

    @validator("*")
    def not_negative(cls, v):
        """Validates that clocks data is not negative."""
        if v < 0.0:
            raise ValueError("Clocks data should be positive.")
        return v


class InteractionsValidator(BaseModel):
    """Class to parse and validate data on object interactions."""

    type: str
    sphere_sphere_adhesion: float
    sphere_sphere_repulsion: float
    sphere_cylinder_adhesion: float
    sphere_cylinder_repulsion: float
    cylinder_cylinder_adhesion: float
    cylinder_cylinder_repulsion: float
    sphere_sphere_smoothness: Optional[int] = None
    sphere_cylinder_smoothness: Optional[int] = None
    cylinder_cylinder_smoothness: Optional[int] = None

    @validator("sphere_sphere_adhesion")
    def not_negative(cls, v):
        """Validates that ineractions data is not negative."""
        if v < 0.0:
            raise ValueError("Physical interactions data should be positive.")
        return v


class ConfigParser:
    """
    Class to read the inputs from the YAML file and validate them.

    Parameters
    ----------
    config_path
        The path to the YAML file to be parsed.
    """

    def __init__(self, config_path: Path) -> None:
        with open(config_path) as file:
            self.cfg = yaml.safe_load(file)

    def get_time_data(self):
        """Returns the simulation time data."""
        return dict(TimerValidator(**self.cfg["time"]))

    def get_domain_data(self):
        """Returns the simulation spatial data."""
        return dict(DomainValidator(**self.cfg["domain"]["boundaries"]))

    def get_2d_status(self):
        """Returns the simulation 2D status."""
        return self.cfg["domain"]["use_2d"]

    def get_drag_coefficient(self):
        """Returns the simulation domain's drag coefficient."""
        return self.cfg["domain"]["drag_coefficient"]

    def get_max_number_of_neurites(self):
        """Returns the maximum number of neurites to be set on each neuron."""
        return self.cfg["neurons"]["max_number_of_neurites"]

    def get_objects_data(self):
        """Returns the data to create a factory instance for cell components."""
        return dict(ObjectValidator(**self.cfg["neurons"]["objects"]))

    def get_clocks_data(self):
        """Returns the data to create the biological clocks for the neurons."""
        return dict(ClocksValidator(**self.cfg["neurons"]["clocks"]))

    def get_interactions_data(self):
        """Returns the data to create object physical interactions."""
        interactions = dict(InteractionsValidator(**self.cfg["interactions"]))
        return {k: v for k, v in interactions.items() if v is not None}
