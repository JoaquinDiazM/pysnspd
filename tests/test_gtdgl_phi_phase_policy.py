"""Smoke tests for selectable OE7 temporal gauge-link policies."""
from __future__ import annotations

import numpy as np
import pytest

from pysnspd.gtdgl.kwt_update import VALID_PHI_PHASE_POLICIES, kwt_delta_update_attempt


class _Defs:
    alpha_kwt_J_inv2 = 0.0
    rho = np.ones(3)
    forcing_J = np.zeros(3, dtype=np.complex128)


class _Material:
    tau0_GL_s = 1.0


def test_phi_phase_policy_names_are_stable():
    assert VALID_PHI_PHASE_POLICIES == {"plus", "none", "minus"}


def test_phi_phase_policy_none_preserves_zero_forcing_state():
    psi = np.array([1.0 + 0.0j, 0.3 + 0.4j, -0.2 + 0.1j])
    phi = np.array([1.0, -2.0, 0.5])
    out, discr = kwt_delta_update_attempt(
        psi_J=psi,
        phi_V=phi,
        defs=_Defs(),
        dt_s=1.0e-16,
        material=_Material(),
        phi_phase_policy="none",
    )
    assert discr >= 0.0
    np.testing.assert_allclose(out, psi)


def test_phi_phase_policy_invalid_raises():
    with pytest.raises(ValueError):
        kwt_delta_update_attempt(
            psi_J=np.ones(2, dtype=np.complex128),
            phi_V=np.zeros(2),
            defs=_Defs(),
            dt_s=1.0e-16,
            material=_Material(),
            phi_phase_policy="bad",
        )
