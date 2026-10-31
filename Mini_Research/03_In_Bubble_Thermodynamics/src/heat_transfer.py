"""Phenomenological conductance with a signed flux, never clipped."""
import numpy as np


def conductance(R, eta, kg):
    if eta < 0 or not np.isfinite(eta):
        raise ValueError('eta must be finite and nonnegative')
    return 4 * np.pi * eta * kg * R  # W/K


def heat_loss(R, T, gas, eta):
    return conductance(R, eta, gas.kg) * (T - gas.Tl)  # W; negative heats gas


def relaxation_time(R, gas, eta, p):
    return gas.capacity(p) / conductance(R, eta, gas.kg) if eta else np.inf
