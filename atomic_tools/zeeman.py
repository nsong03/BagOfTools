"""Zeeman splitting calculations for atomic physics.

This module provides functions to compute the Zeeman energy shift
for hyperfine levels in atoms. It focuses on Rubidium as a common
example in cold-atom experiments but can be used for arbitrary
atoms when the necessary parameters are provided.
"""

from __future__ import annotations

import numpy as np

# Physical constants (SI units)
MU_B = 9.274009994e-24  # Bohr magneton in J/T
HPLANCK = 6.62607015e-34  # Planck constant in J*s

# Electron spin and orbital g-factors
G_S = 2.00231930436256
G_L = 1.0


def g_j(L: float, S: float, J: float, g_s: float = G_S, g_l: float = G_L) -> float:
    """Landé g-factor for electronic angular momentum J.

    Parameters
    ----------
    L, S, J : float
        Orbital, spin and total electronic angular momentum quantum
        numbers.
    g_s, g_l : float
        Spin and orbital g-factors. Defaults are for an electron.

    Returns
    -------
    float
        Landé g-factor g_J.
    """
    jterm = J * (J + 1)
    lterm = L * (L + 1)
    sterm = S * (S + 1)
    return (
        g_l * (jterm + lterm - sterm) / (2 * jterm)
        + g_s * (jterm + sterm - lterm) / (2 * jterm)
    )


def g_f(I: float, J: float, F: float, g_j_val: float, g_i: float) -> float:
    """Hyperfine Landé g-factor g_F.

    Parameters
    ----------
    I, J, F : float
        Nuclear, electronic and total angular momenta.
    g_j_val : float
        Electronic Landé g-factor g_J.
    g_i : float
        Nuclear g-factor.
    Returns
    -------
    float
        Landé g-factor g_F for the hyperfine level.
    """
    fterm = F * (F + 1)
    jterm = J * (J + 1)
    iterm = I * (I + 1)
    return (
        g_j_val * (fterm + jterm - iterm) / (2 * fterm)
        + g_i * (fterm + iterm - jterm) / (2 * fterm)
    )


def zeeman_shift(
    B: np.ndarray | float,
    m_f: np.ndarray | float,
    g_f_val: float,
    mu_b: float = MU_B,
) -> np.ndarray:
    """Zeeman energy shift ΔE for given magnetic field and m_F levels.

    Parameters
    ----------
    B : array_like or float
        Magnetic field magnitude in Tesla.
    m_f : array_like or float
        Magnetic quantum number(s) m_F.
    g_f_val : float
        Landé g-factor g_F for the hyperfine level.
    mu_b : float, optional
        Bohr magneton to use, defaults to MU_B.

    Returns
    -------
    numpy.ndarray
        Energy shift(s) ΔE in Joules.
    """
    B = np.asarray(B, dtype=float)
    m_f = np.asarray(m_f, dtype=float)
    return mu_b * g_f_val * B * m_f


def zeeman_frequency(
    B: np.ndarray | float,
    m_f: np.ndarray | float,
    g_f_val: float,
    mu_b: float = MU_B,
    h: float = HPLANCK,
) -> np.ndarray:
    """Zeeman frequency shift Δν = ΔE / h.

    Parameters
    ----------
    B, m_f, g_f_val : see :func:`zeeman_shift`.
    mu_b : float, optional
        Bohr magneton.
    h : float, optional
        Planck constant.

    Returns
    -------
    numpy.ndarray
        Frequency shift(s) Δν in Hz.
    """
    energy = zeeman_shift(B, m_f, g_f_val, mu_b)
    return energy / h


# Convenience functions for Rubidium ground state

# Nuclear spins
I_RB87 = 3.0 / 2.0
I_RB85 = 5.0 / 2.0
# Nuclear g-factors (CODATA 2018)
G_I_RB87 = -0.0009951414
G_I_RB85 = -0.0002936400


def rb87_ground_gf(F: int) -> float:
    """Return g_F for 87Rb 5S1/2 ground state.

    Parameters
    ----------
    F : int
        Hyperfine level, 1 or 2.

    Returns
    -------
    float
        Landé g-factor g_F.
    """
    J = 0.5
    L = 0.0
    S = 0.5
    gj = g_j(L, S, J)
    return g_f(I_RB87, J, float(F), gj, G_I_RB87)


def rb85_ground_gf(F: int) -> float:
    """Return g_F for 85Rb 5S1/2 ground state (F=2 or 3)."""
    J = 0.5
    L = 0.0
    S = 0.5
    gj = g_j(L, S, J)
    return g_f(I_RB85, J, float(F), gj, G_I_RB85)


