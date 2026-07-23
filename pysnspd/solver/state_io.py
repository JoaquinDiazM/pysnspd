"""NPZ persistence helpers for OE7 stationary gTDGL state/history."""
from __future__ import annotations

from pathlib import Path
import json

import numpy as np

from pysnspd.gtdgl.state import GTDGLStationaryState


def save_stationary_state_npz(state: GTDGLStationaryState, output_path: str | Path) -> Path:
    """Save a relaxed stationary state to NPZ."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    c = state.currents
    np.savez_compressed(
        output,
        psi_real_J=np.real(state.psi_J),
        psi_imag_J=np.imag(state.psi_J),
        phi_V=state.phi_V,
        Te_K=state.Te_K,
        Tph_K=state.Tph_K,
        edge_Q_m_inv=c.edge_Q_m_inv,
        edge_js_us_A_m2=c.edge_js_us_A_m2,
        edge_js_gl_A_m2=c.edge_js_gl_A_m2,
        edge_jn_A_m2=c.edge_jn_A_m2,
        edge_jtot_A_m2=c.edge_jtot_A_m2,
        node_div_js_us_A_m3=c.node_div_js_us_A_m3,
        node_div_js_gl_A_m3=c.node_div_js_gl_A_m3,
        node_div_jtot_A_m3=c.node_div_jtot_A_m3,
        node_js_us_x_A_m2=c.node_js_us_x_A_m2,
        node_js_us_y_A_m2=c.node_js_us_y_A_m2,
        node_jn_x_A_m2=c.node_jn_x_A_m2,
        node_jn_y_A_m2=c.node_jn_y_A_m2,
        node_jtot_x_A_m2=c.node_jtot_x_A_m2,
        node_jtot_y_A_m2=c.node_jtot_y_A_m2,
        edge_pairbreaking_ratio=c.edge_pairbreaking_ratio,
        node_pairbreaking_ratio=c.node_pairbreaking_ratio,
        metadata_json=json.dumps(state.metadata, sort_keys=True),
    )
    return output


def save_relaxation_history_npz(history: dict[str, np.ndarray], output_path: str | Path) -> Path:
    """Save compact relaxation history to NPZ."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    arrays = {key: np.asarray(value) for key, value in history.items()}
    np.savez_compressed(output, **arrays)
    return output

