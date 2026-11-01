"""One optical absorption and one signed Au-water channel per local particle."""
from dataclasses import dataclass, asdict
from functools import cached_property
import numpy as np
from .upstream import optics

@dataclass(frozen=True)
class Particle:
    count: int = 1
    radius_m: float = 20e-9
    gap_m: float = 100e-9
    motion: str = 'capture'
    wavelength_m: float = 532e-9
    G_W_m2K: float = 1e8
    beam_radius_m: float = 100e-6
    absorption_scale: float = 1.

    def __post_init__(self):
        if self.count<0 or int(self.count)!=self.count or min(self.radius_m,self.gap_m,self.beam_radius_m)<=0 or min(self.G_W_m2K,self.absorption_scale)<0:
            raise ValueError('Invalid particle/optical inputs')
        if self.motion not in ('capture','fixed','passive'): raise ValueError('Unknown motion')
        if self.absorption_scale>1:raise ValueError('External absorption allowance must be in [0,1]')

    @cached_property
    def capacity(self):
        return optics.thermal_scales(self.radius_m)['capacity_J_K']

    def optical(self):
        c=optics.optical_cross_sections(self.radius_m,self.wavelength_m,n_medium=1.33)
        if self.count*c.extinction_m2>np.pi*self.beam_radius_m**2:
            raise ValueError('Local extinction exceeds beam allowance')
        return c

    def metadata(self):
        c=self.optical()
        return {**asdict(self),'sigma_abs_m2':float(c.absorption_m2),
                'sigma_ext_m2':float(c.extinction_m2),'capacity_J_K':self.capacity,
                'N_sigma_ext_over_beam_area':float(self.count*c.extinction_m2/(np.pi*self.beam_radius_m**2))}
