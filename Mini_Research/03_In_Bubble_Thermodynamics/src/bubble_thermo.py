"""Ideal-gas subsystem; vapor is excluded from n, pressure and internal energy."""
from dataclasses import dataclass
from pathlib import Path
import importlib.util
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
MR2 = ROOT.parent / '02_Compressible_Bubble_Dynamics'
MECHANICAL_COMMIT = 'f1f3e812d77fd7fa058ce111ff484b1c8658f937'
MR1_COMMIT = 'd59bc0d0119b514152f763e2537a0477c7c63663'
spec = importlib.util.spec_from_file_location('mr3_accepted_mechanics', MR2 / 'src/model.py')
mechanics = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mechanics
spec.loader.exec_module(mechanics)
Parameters = mechanics.Parameters
RU = 8.31446261815324  # J mol^-1 K^-1; exact SI molar gas constant


@dataclass(frozen=True)
class Gas:
    gamma: float = 1.4
    T0: float = 293.15
    Tl: float = 293.15
    kg: float = 0.0264  # NIST Huber table: AIR, 300 K, 100 kPa
    molar_mass: float = 0.02897  # kg/mol, idealized dry-air bookkeeping
    name: str = 'air_like'

    def __post_init__(self):
        if not (self.gamma > 1 and min(self.T0, self.Tl, self.kg, self.molar_mass) > 0
                and np.all(np.isfinite([self.gamma, self.T0, self.Tl, self.kg, self.molar_mass]))):
            raise ValueError('Finite positive gas parameters and gamma > 1 required')

    @property
    def Cv(self):
        return RU / (self.gamma - 1)

    def n(self, p=Parameters()):
        return p.pg0 * volume(p.R0) / (RU * self.T0)

    def capacity(self, p=Parameters()):
        return self.n(p) * self.Cv

    def E0(self, p=Parameters()):
        return self.capacity(p) * self.T0


def volume(R):
    return 4 * np.pi * np.asarray(R)**3 / 3


def adiabatic_temperature(R, gas=Gas(), p=Parameters()):
    return gas.T0 * (p.R0 / np.asarray(R))**(3 * (gas.gamma - 1))


def adiabatic_pressure(R, gas=Gas(), p=Parameters()):
    return p.pg0 * (p.R0 / np.asarray(R))**(3 * gas.gamma)


def gas_pressure(R, T, gas=Gas(), p=Parameters()):
    return gas.n(p) * RU * np.asarray(T) / volume(R)


def gas_density(R, gas=Gas(), p=Parameters()):
    return gas.n(p) * gas.molar_mass / volume(R)


def temperature_rate(R, U, T, qnet, gas=Gas(), p=Parameters()):
    return -3 * (gas.gamma - 1) * U * T / R + qnet / gas.capacity(p)


def pressure_rate(R, U, T, qnet, gas=Gas(), p=Parameters()):
    pg = gas_pressure(R, T, gas, p)
    return -3 * gas.gamma * pg * U / R + (gas.gamma - 1) * qnet / volume(R)
