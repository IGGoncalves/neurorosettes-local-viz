import yaml
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, validator


class TimerValidator(BaseModel):
    total_time: float
    step: float

    @validator("total_time", "step")
    def not_negative(cls, v):
        if v < 0.0:
            raise ValueError("Time data should be positive.")
        return v


class DomainValidator(BaseModel):
    min: float
    max: float
    step: float

    @validator("step")
    def not_negative(cls, v):
        if v < 0.0:
            raise ValueError("Domain step should be positive.")
        return v


class ObjectValidator(BaseModel):
    cell_radius: float
    cell_interaction_factor: float
    neurite_radius: float
    neurite_interaction_factor: float
    neurite_spring_constant: float
    neurite_default_length: float

    @validator("cell_radius")
    def not_negative(cls, v):
        if v < 0.0:
            raise ValueError("Physical object data should be positive.")
        return v

class ClocksValidator(BaseModel):
    proliferation_rate: float
    death_rate: float
    differentiation_rate: float

    @validator("*")
    def not_negative(cls, v):
        if v < 0.0:
            raise ValueError("Physical object data should be positive.")
        return v


class ClocksValidator(BaseModel):
    proliferation_rate: float
    death_rate: float
    differentiation_rate: float

class InteractionsValidator(BaseModel):
    type: str
    sphere_sphere_adhesion: float
    sphere_sphere_repulsion: float
    sphere_sphere_smoothness: Optional[int]
    sphere_cylinder_adhesion: float
    sphere_cylinder_repulsion: float
    sphere_cylinder_smoothness: Optional[int]
    cylinder_cylinder_adhesion: float
    cylinder_cylinder_repulsion: float
    cylinder_cylinder_smoothness: Optional[int]

    @validator("sphere_sphere_adhesion")
    def not_negative(cls, v):
        if v < 0.0:
            raise ValueError("Physical object data should be positive.")
        return v


class ConfigParser:
    def __init__(self, config_path: Path) -> None:
        with open(config_path) as file:
            self.cfg = yaml.safe_load(file)

    def get_time_data(self):
        return dict(TimerValidator(**self.cfg["time"]))

    def get_domain_data(self):
        return dict(DomainValidator(**self.cfg["domain"]["boundaries"]))

    def get_2d_status(self):
        return self.cfg["domain"]["use_2d"]

    def get_drag_coefficient(self):
        return self.cfg["domain"]["drag_coefficient"]

    def get_max_number_of_neurites(self):
        return self.cfg["neurons"]["max_number_of_neurites"]

    def get_objects_data(self):
        return dict(ObjectValidator(**self.cfg["neurons"]["objects"]))

    def get_clocks_data(self):
        return dict(ClocksValidator(**self.cfg["neurons"]["clocks"]))

    def get_interactions_data(self):
        return dict(InteractionsValidator(**self.cfg["interactions"]))
