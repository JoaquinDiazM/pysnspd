"""A-posteriori energetic projection diagnostics for persisted photon runs.

This module deliberately does not advance the coupled solver.  It streams the
large compressed snapshot archive, projects the already-persisted trajectory on
the Simon electronic-energy catalogue, and estimates the spectral-storage
terms omitted by the present temperature update.
"""

from __future__ import annotations

from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
import time
from typing import Any, Callable, Iterator, Mapping, Sequence
from zipfile import ZipFile

import numpy as np

from pysnspd.analysis.photon_snapshots import _SnapshotProjector
from pysnspd.thermal.evolution import _interp_nd

MEV_J = 1.602176634e-22


@dataclass(frozen=True)
class EnergyProjectionCatalog:
    """Minimal view of the PRE catalogue needed by the D3 diagnostic."""

    Te_values_K: np.ndarray
    Tph_values_K: np.ndarray
    delta_values_J: np.ndarray
    q_values_m_inv: np.ndarray
    u_e_J_m3: np.ndarray
    C_e_J_m3_K: np.ndarray
    P_total_W_m3: np.ndarray
    kappa_s_W_m_K: np.ndarray
    u_ph_J_m3: np.ndarray
    P_esc_W_m3: np.ndarray
    N0_J_m3: float
    delta0_J: float
    bath_K: float

    @classmethod
    def load(cls, path: str | Path) -> "EnergyProjectionCatalog":
        with np.load(Path(path), allow_pickle=True) as data:
            metadata = dict(np.asarray(data["metadata"], dtype=object).item())
            return cls(
                Te_values_K=np.asarray(data["Te_values_K"], dtype=float),
                Tph_values_K=np.asarray(data["Tph_values_K"], dtype=float),
                delta_values_J=np.asarray(data["delta_values_J"], dtype=float),
                q_values_m_inv=np.asarray(data["q_values_m_inv"], dtype=float),
                u_e_J_m3=np.asarray(data["u_e_J_m3"], dtype=float),
                C_e_J_m3_K=np.asarray(data["C_e_J_m3_K"], dtype=float),
                P_total_W_m3=np.asarray(data["P_total_W_m3"], dtype=float),
                kappa_s_W_m_K=np.asarray(data["kappa_s_W_m_K"], dtype=float),
                u_ph_J_m3=np.asarray(data["u_ph_J_m3"], dtype=float),
                P_esc_W_m3=np.asarray(data["P_esc_W_m3"], dtype=float),
                N0_J_m3=float(metadata["N0_J_m3"]),
                delta0_J=float(metadata["delta0_J"]),
                bath_K=float(metadata.get("T_bath_K", np.nan)),
            )

    def electronic_energy(
        self,
        Te_K: np.ndarray,
        delta_J: np.ndarray,
        q_abs_m_inv: np.ndarray,
    ) -> np.ndarray:
        coords = self.clip_state(Te_K, None, delta_J, q_abs_m_inv)
        return np.asarray(
            _interp_nd(
                self.u_e_J_m3,
                (self.Te_values_K, self.delta_values_J, self.q_values_m_inv),
                (coords[0], coords[2], coords[3]),
            ),
            dtype=float,
        )

    def evaluate(
        self,
        *,
        Te_K: np.ndarray,
        Tph_K: np.ndarray,
        delta_J: np.ndarray,
        q_abs_m_inv: np.ndarray,
    ) -> dict[str, np.ndarray]:
        Te, Tph, delta, q_abs = self.clip_state(Te_K, Tph_K, delta_J, q_abs_m_inv)
        axes3 = (self.Te_values_K, self.delta_values_J, self.q_values_m_inv)
        coords3 = (Te, delta, q_abs)
        axes4 = (
            self.Te_values_K,
            self.Tph_values_K,
            self.delta_values_J,
            self.q_values_m_inv,
        )
        coords4 = (Te, Tph, delta, q_abs)
        axes2 = (self.Te_values_K, self.delta_values_J)
        coords2 = (Te, delta)
        return {
            "u_e_J_m3": np.asarray(_interp_nd(self.u_e_J_m3, axes3, coords3), dtype=float),
            "C_e_J_m3_K": np.asarray(_interp_nd(self.C_e_J_m3_K, axes3, coords3), dtype=float),
            "P_ep_W_m3": np.asarray(_interp_nd(self.P_total_W_m3, axes4, coords4), dtype=float),
            "kappa_W_m_K": np.asarray(_interp_nd(self.kappa_s_W_m_K, axes2, coords2), dtype=float),
            "u_ph_J_m3": np.asarray(np.interp(Tph, self.Tph_values_K, self.u_ph_J_m3), dtype=float),
            "P_esc_W_m3": np.asarray(np.interp(Tph, self.Tph_values_K, self.P_esc_W_m3), dtype=float),
        }

    def clip_state(
        self,
        Te_K: np.ndarray,
        Tph_K: np.ndarray | None,
        delta_J: np.ndarray,
        q_abs_m_inv: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray | None, np.ndarray, np.ndarray]:
        Te = np.clip(np.asarray(Te_K, dtype=float), self.Te_values_K[0], self.Te_values_K[-1])
        Tph = None
        if Tph_K is not None:
            Tph = np.clip(
                np.asarray(Tph_K, dtype=float),
                self.Tph_values_K[0],
                self.Tph_values_K[-1],
            )
        delta = np.clip(
            np.asarray(delta_J, dtype=float),
            self.delta_values_J[0],
            self.delta_values_J[-1],
        )
        q_abs = np.clip(
            np.asarray(q_abs_m_inv, dtype=float),
            self.q_values_m_inv[0],
            self.q_values_m_inv[-1],
        )
        return Te, Tph, delta, q_abs

    def clipping_mask(
        self,
        Te_K: np.ndarray,
        Tph_K: np.ndarray,
        delta_J: np.ndarray,
        q_abs_m_inv: np.ndarray,
    ) -> np.ndarray:
        components = self.clipping_components(Te_K, Tph_K, delta_J, q_abs_m_inv)
        return np.logical_or.reduce(tuple(components.values()))

    def clipping_components(
        self,
        Te_K: np.ndarray,
        Tph_K: np.ndarray,
        delta_J: np.ndarray,
        q_abs_m_inv: np.ndarray,
    ) -> dict[str, np.ndarray]:
        return {
            "Te": (np.asarray(Te_K) < self.Te_values_K[0])
            | (np.asarray(Te_K) > self.Te_values_K[-1]),
            "Tph": (np.asarray(Tph_K) < self.Tph_values_K[0])
            | (np.asarray(Tph_K) > self.Tph_values_K[-1]),
            "delta": (np.asarray(delta_J) < self.delta_values_J[0])
            | (np.asarray(delta_J) > self.delta_values_J[-1]),
            "q": (np.asarray(q_abs_m_inv) < self.q_values_m_inv[0])
            | (np.asarray(q_abs_m_inv) > self.q_values_m_inv[-1]),
        }

    def condensation_energy(self, delta_J: np.ndarray) -> np.ndarray:
        delta = np.maximum(np.asarray(delta_J, dtype=float), 0.0)
        result = np.zeros_like(delta)
        positive = delta > 0.0
        result[positive] = -self.N0_J_m3 * delta[positive] ** 2 * (
            0.5 + np.log(self.delta0_J / np.maximum(delta[positive], 1.0e-300))
        )
        return result


@dataclass(frozen=True)
class SupercurrentProjectionCatalog:
    """Strict PRE supercurrent table used to audit out-of-axis momentum states."""

    Te_axis_K: np.ndarray
    delta_axis_J: np.ndarray
    q_axis_m_inv: np.ndarray
    js_A_m2: np.ndarray

    @classmethod
    def load(cls, path: str | Path) -> "SupercurrentProjectionCatalog":
        with np.load(Path(path), allow_pickle=True) as data:
            Te_axis = _first_npz_array(data, ("Te_axis_K", "Te_values_K", "T_axis_K"))
            delta_axis = _first_npz_array(
                data,
                ("delta_axis_J", "delta_values_J", "Delta_axis_J"),
            )
            q_axis = _first_npz_array(data, ("q_axis_m_inv", "q_values_m_inv"))
            table = _first_npz_array(data, ("js_A_m2", "j_s_A_m2", "js_T_delta_q_A_m2"))
        expected = (Te_axis.size, delta_axis.size, q_axis.size)
        if table.shape != expected:
            raise ValueError(
                "strict supercurrent table must have layout js_A_m2[Te,delta,q] "
                f"with shape {expected}; received {table.shape}"
            )
        if any(axis.size == 0 or np.any(~np.isfinite(axis)) for axis in (Te_axis, delta_axis, q_axis)):
            raise ValueError("strict supercurrent axes must be finite and nonempty")
        if any(np.any(np.diff(axis) < 0.0) for axis in (Te_axis, delta_axis, q_axis)):
            raise ValueError("strict supercurrent axes must be monotonically increasing")
        if np.any(~np.isfinite(table)):
            raise ValueError("strict supercurrent table contains non-finite values")
        return cls(
            Te_axis_K=Te_axis,
            delta_axis_J=delta_axis,
            q_axis_m_inv=q_axis,
            js_A_m2=table,
        )

    @property
    def q_max_m_inv(self) -> float:
        return float(self.q_axis_m_inv[-1])

    @property
    def js_reference_A_m2(self) -> float:
        return float(np.nanmax(np.abs(self.js_A_m2)))

    def evaluate_abs(
        self,
        Te_K: np.ndarray,
        delta_J: np.ndarray,
        q_abs_m_inv: np.ndarray,
    ) -> np.ndarray:
        Te = np.clip(np.asarray(Te_K, dtype=float), self.Te_axis_K[0], self.Te_axis_K[-1])
        delta = np.clip(
            np.asarray(delta_J, dtype=float),
            self.delta_axis_J[0],
            self.delta_axis_J[-1],
        )
        q_abs = np.clip(
            np.asarray(q_abs_m_inv, dtype=float),
            self.q_axis_m_inv[0],
            self.q_axis_m_inv[-1],
        )
        return np.abs(
            np.asarray(
                _interp_nd(
                    self.js_A_m2,
                    (self.Te_axis_K, self.delta_axis_J, self.q_axis_m_inv),
                    (Te, delta, q_abs),
                ),
                dtype=float,
            )
        )


class CompressedNpzRowStream:
    """Read selected 2D ``.npy`` members of a compressed NPZ by row chunks.

    ``numpy.load`` inflates each requested member completely.  The photon file
    contains several 1.15-GiB members, so keeping independent sequential ZIP
    streams bounds memory by ``chunk_rows * n_nodes`` instead.
    """

    def __init__(self, path: str | Path, keys: Sequence[str]):
        self.path = Path(path)
        self.keys = tuple(str(key) for key in keys)
        if not self.keys:
            raise ValueError("At least one NPZ member is required.")

    def iter_chunks(
        self,
        *,
        chunk_rows: int,
        stop_row: int | None = None,
    ) -> Iterator[tuple[int, dict[str, np.ndarray]]]:
        stride = max(1, int(chunk_rows))
        with ExitStack() as stack:
            streams: dict[str, Any] = {}
            specs: dict[str, tuple[np.dtype, tuple[int, ...]]] = {}
            for key in self.keys:
                archive = stack.enter_context(ZipFile(self.path))
                stream = stack.enter_context(archive.open(f"{key}.npy", "r"))
                version = np.lib.format.read_magic(stream)
                if version == (1, 0):
                    shape, fortran, dtype = np.lib.format.read_array_header_1_0(stream)
                elif version == (2, 0):
                    shape, fortran, dtype = np.lib.format.read_array_header_2_0(stream)
                else:
                    shape, fortran, dtype = np.lib.format._read_array_header(stream, version)
                if fortran or len(shape) != 2 or np.dtype(dtype).hasobject:
                    raise ValueError(
                        f"{key} must be a C-order, non-object 2D array; got shape={shape}."
                    )
                streams[key] = stream
                specs[key] = (np.dtype(dtype), tuple(int(value) for value in shape))

            shapes = {shape for _, shape in specs.values()}
            if len(shapes) != 1:
                raise ValueError(f"Streamed NPZ members do not share a shape: {specs}.")
            n_rows, n_columns = next(iter(shapes))
            stop = n_rows if stop_row is None else min(n_rows, max(0, int(stop_row)))
            offset = 0
            while offset < stop:
                count = min(stride, stop - offset)
                chunk: dict[str, np.ndarray] = {}
                for key, stream in streams.items():
                    dtype, _ = specs[key]
                    nbytes = count * n_columns * dtype.itemsize
                    raw = _read_exact(stream, nbytes)
                    chunk[key] = np.frombuffer(raw, dtype=dtype).reshape(count, n_columns)
                yield offset, chunk
                offset += count


def extract_energy_projection_diagnostics(
    *,
    snapshots_npz: str | Path,
    power_table_npz: str | Path,
    usadel_current_npz: str | Path,
    history: Mapping[str, Any],
    nodes_m: np.ndarray,
    triangles: np.ndarray,
    ops: Any,
    sigma_n_S_m: float,
    thickness_m: float,
    Tc_K: float,
    xi_m: float,
    window_m: float = 100.0e-9,
    chunk_rows: int = 32,
    snapshot_stride: int = 1,
    max_snapshots: int | None = None,
    progress_callback: Callable[[int, int, float], None] | None = None,
) -> dict[str, np.ndarray]:
    """Stream a photon trajectory and return a compact D3 plotting dataset."""

    nodes = np.asarray(nodes_m, dtype=float)
    triangles_array = np.asarray(triangles, dtype=np.int64)
    if nodes.ndim != 2 or nodes.shape[1] < 2:
        raise ValueError("nodes_m must have shape (n_nodes, >=2).")
    if float(thickness_m) <= 0.0 or float(window_m) <= 0.0:
        raise ValueError("thickness_m and window_m must be positive.")
    if float(xi_m) <= 0.0 or float(Tc_K) <= 0.0:
        raise ValueError("xi_m and Tc_K must be positive.")

    snapshot_times_ps = _load_snapshot_times(snapshots_npz)
    total_snapshots = int(snapshot_times_ps.size)
    stop = total_snapshots if max_snapshots is None else min(total_snapshots, int(max_snapshots))
    if stop < 2:
        raise ValueError("At least two persisted snapshots are required.")
    raw_times = snapshot_times_ps[:stop]
    sample_stride = max(1, int(snapshot_stride))
    candidate_indices = np.arange(0, stop, sample_stride, dtype=np.int64)
    if candidate_indices[-1] != stop - 1:
        candidate_indices = np.append(candidate_indices, stop - 1)
    sampled_indices = _strictly_increasing_sample_indices(
        raw_times,
        candidate_indices,
    )
    dropped_duplicate_times = int(candidate_indices.size - sampled_indices.size)
    sampled_lookup = set(int(value) for value in sampled_indices)

    history_t = _history_time_ps(history)
    history_v = np.asarray(history.get("V_tdgl_center_V", []), dtype=float).reshape(-1)
    hcount = min(history_t.size, history_v.size)
    valid = np.isfinite(history_t[:hcount]) & np.isfinite(history_v[:hcount]) & (history_t[:hcount] <= raw_times[-1] + 1.0e-9)
    if np.any(valid):
        candidates = np.flatnonzero(valid)
        vmax_history_index = int(candidates[np.argmax(history_v[:hcount][valid])])
        vmax_time_ps = float(history_t[vmax_history_index])
    else:
        vmax_time_ps = float(raw_times[0])
    final_time_ps = float(raw_times[-1])
    midpoint_time_ps = 0.5 * (vmax_time_ps + final_time_ps)
    requested_map_times = np.asarray(
        [vmax_time_ps, midpoint_time_ps, final_time_ps], dtype=float
    )
    target_indices = np.asarray(
        [sampled_indices[np.argmin(np.abs(raw_times[sampled_indices] - target))] for target in requested_map_times],
        dtype=np.int64,
    )
    target_indices = np.unique(target_indices)
    target_set = set(int(value) for value in target_indices if int(value) > 0)

    x = nodes[:, 0]
    center_x = 0.5 * (float(np.nanmin(x)) + float(np.nanmax(x)))
    central_mask = np.abs(x - center_x) <= 0.5 * float(window_m)
    if not np.any(central_mask):
        raise ValueError("The requested central window contains no mesh nodes.")
    center_indices = np.flatnonzero(central_mask)
    node_area = np.asarray(ops.node_area_m2, dtype=float).reshape(-1)
    if node_area.size != nodes.shape[0]:
        raise ValueError("Finite-volume node areas do not match the mesh.")
    volume_weights = node_area[center_indices] * float(thickness_m)
    volume_total = float(np.sum(volume_weights))

    catalog = EnergyProjectionCatalog.load(power_table_npz)
    current_catalog = SupercurrentProjectionCatalog.load(usadel_current_npz)
    projector = _SnapshotProjector(ops)
    bath_K = float(catalog.bath_K)
    if not np.isfinite(bath_K) or bath_K <= 0.0:
        bath_K = 0.9

    temporal_lists: dict[str, list[np.ndarray]] = {
        key: [] for key in _TEMPORAL_KEYS
    }
    selected_records: dict[int, dict[str, np.ndarray]] = {}
    initial_maps: dict[str, np.ndarray] | None = None
    carry: dict[str, np.ndarray] | None = None
    carry_indices = np.array([], dtype=np.int64)
    started = time.monotonic()

    stream = CompressedNpzRowStream(
        snapshots_npz,
        ("psi_real_snapshot_J", "psi_imag_snapshot_J", "phi_snapshot_V", "Te_snapshot_K", "Tph_snapshot_K"),
    )
    for offset, raw_chunk in stream.iter_chunks(
        chunk_rows=max(1, int(chunk_rows)) * sample_stride,
        stop_row=stop,
    ):
        raw_indices = np.arange(offset, offset + next(iter(raw_chunk.values())).shape[0], dtype=np.int64)
        keep = np.asarray([int(index) in sampled_lookup for index in raw_indices], dtype=bool)
        if not np.any(keep):
            if progress_callback is not None:
                progress_callback(int(raw_indices[-1]) + 1, stop, time.monotonic() - started)
            continue
        chunk = {key: np.asarray(value[keep], dtype=float) for key, value in raw_chunk.items()}
        indices = raw_indices[keep]
        if carry is not None:
            chunk = {
                key: np.concatenate([carry[key], value], axis=0)
                for key, value in chunk.items()
            }
            indices = np.concatenate([carry_indices, indices])

        derived = _derive_states(
            raw=chunk,
            catalog=catalog,
            current_catalog=current_catalog,
            projector=projector,
            center_indices=center_indices,
            central_mask=central_mask,
            sigma_n_S_m=float(sigma_n_S_m),
            bath_K=bath_K,
        )
        if initial_maps is None:
            initial_maps = _full_state_maps(
                position=0,
                raw=chunk,
                derived=derived,
                catalog=catalog,
                xi_m=float(xi_m),
            )

        interval = _interval_fields(
            raw=chunk,
            derived=derived,
            indices=indices,
            times_ps=raw_times,
            catalog=catalog,
            center_indices=center_indices,
        )
        if indices.size >= 2:
            _append_temporal_metrics(
                temporal_lists,
                interval=interval,
                derived=derived,
                center_indices=center_indices,
                volume_weights=volume_weights,
                volume_total=volume_total,
                xi_m=float(xi_m),
                delta0_J=float(catalog.delta0_J),
            )
            for local_end in range(1, indices.size):
                global_end = int(indices[local_end])
                if global_end in target_set:
                    selected_records[global_end] = _selected_interval_maps(
                        local_end=local_end,
                        raw=chunk,
                        derived=derived,
                        interval=interval,
                        catalog=catalog,
                        projector=projector,
                        central_mask=central_mask,
                        sigma_n_S_m=float(sigma_n_S_m),
                        bath_K=bath_K,
                        xi_m=float(xi_m),
                    )

        carry = {key: np.asarray(value[-1:], dtype=float).copy() for key, value in chunk.items()}
        carry_indices = np.asarray(indices[-1:], dtype=np.int64)
        if progress_callback is not None:
            progress_callback(int(raw_indices[-1]) + 1, stop, time.monotonic() - started)

    if initial_maps is None or not temporal_lists["time_ps"]:
        raise RuntimeError("The streamed trajectory did not produce diagnostic intervals.")
    missing_targets = target_set.difference(selected_records)
    if missing_targets:
        raise RuntimeError(f"Selected diagnostic maps were not extracted: {sorted(missing_targets)}")

    selected_sorted = sorted(selected_records)
    result: dict[str, np.ndarray] = {
        "nodes_x_nm": 1.0e9 * nodes[:, 0],
        "nodes_y_nm": 1.0e9 * nodes[:, 1],
        "triangles": triangles_array,
        "node_area_m2": node_area,
        "central_mask": central_mask,
        "selected_snapshot_indices": np.asarray(selected_sorted, dtype=np.int64),
        "selected_times_ps": raw_times[np.asarray(selected_sorted, dtype=np.int64)],
        "requested_selected_times_ps": requested_map_times,
        "photon_time_ps": np.asarray([_photon_time_ps(history)], dtype=float),
        "vmax_time_ps": np.asarray([vmax_time_ps], dtype=float),
        "final_time_ps": np.asarray([final_time_ps], dtype=float),
        "window_nm": np.asarray([1.0e9 * float(window_m)], dtype=float),
        "xi_m": np.asarray([float(xi_m)], dtype=float),
        "Tc_K": np.asarray([float(Tc_K)], dtype=float),
        "delta0_J": np.asarray([float(catalog.delta0_J)], dtype=float),
        "strict_q_max_m_inv": np.asarray([current_catalog.q_max_m_inv], dtype=float),
        "strict_js_reference_A_m2": np.asarray(
            [current_catalog.js_reference_A_m2], dtype=float
        ),
        "thickness_m": np.asarray([float(thickness_m)], dtype=float),
        "total_snapshot_count": np.asarray([total_snapshots], dtype=np.int64),
        "processed_snapshot_count": np.asarray([stop], dtype=np.int64),
        "snapshot_stride": np.asarray([sample_stride], dtype=np.int64),
        "dropped_duplicate_time_count": np.asarray(
            [dropped_duplicate_times], dtype=np.int64
        ),
        "truncated": np.asarray([stop < total_snapshots], dtype=bool),
    }
    for key, values in temporal_lists.items():
        result[key] = np.concatenate(values, axis=0)
    for key in initial_maps:
        result[f"initial_{key}"] = np.asarray(initial_maps[key], dtype=float)
    for key in next(iter(selected_records.values())):
        result[f"selected_{key}"] = np.stack(
            [selected_records[index][key] for index in selected_sorted], axis=0
        )

    temporal_time = result["time_ps"]
    if history_t.size and history_v.size:
        count = min(history_t.size, history_v.size)
        finite = np.isfinite(history_t[:count]) & np.isfinite(history_v[:count])
        result["V_tdgl_V"] = np.interp(
            temporal_time,
            history_t[:count][finite],
            history_v[:count][finite],
        )
    else:
        result["V_tdgl_V"] = np.full_like(temporal_time, np.nan)
    return result


_TEMPORAL_KEYS = (
    "time_ps",
    "dt_ps",
    "p99_abs_P_spec_W_m3",
    "p99_abs_P_delta_W_m3",
    "p99_abs_P_delta_cond_W_m3",
    "p99_abs_P_q_W_m3",
    "p99_abs_P_path_W_m3",
    "p99_abs_Q_ret_W_m3",
    "p99_abs_residual_W_m3",
    "integrated_P_spec_W",
    "integrated_P_delta_W",
    "integrated_P_delta_cond_W",
    "integrated_P_q_W",
    "integrated_P_path_W",
    "integrated_P_J_W",
    "integrated_minus_P_ep_W",
    "integrated_P_diff_W",
    "integrated_Q_ret_W",
    "integrated_residual_W",
    "integrated_P_esc_W",
    "mean_u_e_J_m3",
    "mean_u_cond_J_m3",
    "mean_u_qp_J_m3",
    "mean_u_ph_J_m3",
    "min_delta_over_delta0",
    "max_Te_K",
    "max_Tph_K",
    "max_q_xi",
    "catalog_clipped_fraction",
    "catalog_Te_clipped_fraction",
    "catalog_Tph_clipped_fraction",
    "catalog_delta_clipped_fraction",
    "catalog_q_clipped_fraction",
    "strict_q_clipped_fraction",
    "strict_q_clipped_js_median_A_m2",
    "strict_q_clipped_js_p95_A_m2",
    "strict_q_clipped_js_max_A_m2",
    "strict_q_clipped_js_p95_over_catalog_max",
    "strict_q_clipped_js_near_zero_fraction",
)


def _derive_states(
    *,
    raw: Mapping[str, np.ndarray],
    catalog: EnergyProjectionCatalog,
    current_catalog: SupercurrentProjectionCatalog,
    projector: _SnapshotProjector,
    center_indices: np.ndarray,
    central_mask: np.ndarray,
    sigma_n_S_m: float,
    bath_K: float,
) -> dict[str, np.ndarray]:
    psi = np.asarray(raw["psi_real_snapshot_J"], dtype=float) + 1j * np.asarray(
        raw["psi_imag_snapshot_J"], dtype=float
    )
    edge_phase = np.angle(psi[:, projector.edge_j] * np.conjugate(psi[:, projector.edge_i]))
    q_full = projector.edge_absolute_to_node(
        np.abs(edge_phase) / np.maximum(projector.edge_length_m[None, :], 1.0e-300)
    )
    delta_full = np.abs(psi)
    Te_full = np.asarray(raw["Te_snapshot_K"], dtype=float)
    Tph_full = np.asarray(raw["Tph_snapshot_K"], dtype=float)
    center_lookup = catalog.evaluate(
        Te_K=Te_full[:, center_indices],
        Tph_K=Tph_full[:, center_indices],
        delta_J=delta_full[:, center_indices],
        q_abs_m_inv=q_full[:, center_indices],
    )
    kappa_full = np.zeros_like(Te_full)
    kappa_full[:, center_indices] = center_lookup["kappa_W_m_K"]
    P_diff_full = projector.diffusion_power_density(
        Te_full,
        kappa_full,
        active_mask=central_mask,
        bath_K=float(bath_K),
    )
    P_J_full = projector.joule_power_density(
        np.asarray(raw["phi_snapshot_V"], dtype=float),
        sigma_n_S_m=float(sigma_n_S_m),
    )
    u_cond_center = catalog.condensation_energy(delta_full[:, center_indices])
    clip_components = catalog.clipping_components(
        Te_full[:, center_indices],
        Tph_full[:, center_indices],
        delta_full[:, center_indices],
        q_full[:, center_indices],
    )
    strict_q_clip_center = q_full[:, center_indices] > current_catalog.q_max_m_inv
    js_abs_center = current_catalog.evaluate_abs(
        Te_full[:, center_indices],
        delta_full[:, center_indices],
        q_full[:, center_indices],
    )
    return {
        "q_full": q_full,
        "delta_full": delta_full,
        "q_center": q_full[:, center_indices],
        "delta_center": delta_full[:, center_indices],
        "Te_center": Te_full[:, center_indices],
        "Tph_center": Tph_full[:, center_indices],
        "u_e_center": center_lookup["u_e_J_m3"],
        "C_e_center": center_lookup["C_e_J_m3_K"],
        "P_ep_center": center_lookup["P_ep_W_m3"],
        "u_ph_center": center_lookup["u_ph_J_m3"],
        "P_esc_center": center_lookup["P_esc_W_m3"],
        "P_diff_center": P_diff_full[:, center_indices],
        "P_J_center": P_J_full[:, center_indices],
        "u_cond_center": u_cond_center,
        "strict_q_clip_center": strict_q_clip_center,
        "js_abs_center_A_m2": js_abs_center,
        "js_reference_A_m2": np.asarray(
            [current_catalog.js_reference_A_m2], dtype=float
        ),
        "clip_center": catalog.clipping_mask(
            Te_full[:, center_indices],
            Tph_full[:, center_indices],
            delta_full[:, center_indices],
            q_full[:, center_indices],
        ),
        **{
            f"clip_{key}_center": np.asarray(values, dtype=bool)
            for key, values in clip_components.items()
        },
    }


def _interval_fields(
    *,
    raw: Mapping[str, np.ndarray],
    derived: Mapping[str, np.ndarray],
    indices: np.ndarray,
    times_ps: np.ndarray,
    catalog: EnergyProjectionCatalog,
    center_indices: np.ndarray,
) -> dict[str, np.ndarray]:
    if indices.size < 2:
        return {}
    dt_s = np.diff(times_ps[indices]) * 1.0e-12
    if np.any(~np.isfinite(dt_s)) or np.any(dt_s <= 0.0):
        raise ValueError("Sampled snapshot times must be finite and strictly increasing.")
    Te = np.asarray(raw["Te_snapshot_K"], dtype=float)[:, center_indices]
    delta = np.asarray(derived["delta_full"], dtype=float)[:, center_indices]
    q_abs = np.asarray(derived["q_full"], dtype=float)[:, center_indices]
    u0 = np.asarray(derived["u_e_center"], dtype=float)[:-1]
    u1 = np.asarray(derived["u_e_center"], dtype=float)[1:]
    u_delta = catalog.electronic_energy(Te[:-1], delta[1:], q_abs[:-1])
    u_new_spectrum = catalog.electronic_energy(Te[:-1], delta[1:], q_abs[1:])
    u_q_first = catalog.electronic_energy(Te[:-1], delta[:-1], q_abs[1:])
    scale = dt_s[:, None]
    P_delta = (u_delta - u0) / scale
    P_q = (u_new_spectrum - u_delta) / scale
    P_spec = (u_new_spectrum - u0) / scale
    P_delta_after_q = (u_new_spectrum - u_q_first) / scale
    P_path = P_delta - P_delta_after_q
    cond = np.asarray(derived["u_cond_center"], dtype=float)
    P_delta_cond = (cond[1:] - cond[:-1]) / scale
    P_thermal = (u1 - u_new_spectrum) / scale
    P_CdT = (
        0.5 * (
            np.asarray(derived["C_e_center"], dtype=float)[:-1]
            + np.asarray(derived["C_e_center"], dtype=float)[1:]
        )
        * (Te[1:] - Te[:-1])
        / scale
    )
    P_J = 0.5 * (
        np.asarray(derived["P_J_center"], dtype=float)[:-1]
        + np.asarray(derived["P_J_center"], dtype=float)[1:]
    )
    P_ep = 0.5 * (
        np.asarray(derived["P_ep_center"], dtype=float)[:-1]
        + np.asarray(derived["P_ep_center"], dtype=float)[1:]
    )
    P_diff = 0.5 * (
        np.asarray(derived["P_diff_center"], dtype=float)[:-1]
        + np.asarray(derived["P_diff_center"], dtype=float)[1:]
    )
    P_esc = 0.5 * (
        np.asarray(derived["P_esc_center"], dtype=float)[:-1]
        + np.asarray(derived["P_esc_center"], dtype=float)[1:]
    )
    Q_ret = P_J + P_diff - P_ep
    du_dt = (u1 - u0) / scale
    return {
        "time_ps": 0.5 * (times_ps[indices[:-1]] + times_ps[indices[1:]]),
        "dt_ps": 1.0e12 * dt_s,
        "P_delta": P_delta,
        "P_delta_cond": P_delta_cond,
        "P_delta_qp": P_delta - P_delta_cond,
        "P_q": P_q,
        "P_spec": P_spec,
        "P_path": P_path,
        "P_thermal": P_thermal,
        "P_CdT": P_CdT,
        "P_T_nonlinear": P_thermal - P_CdT,
        "P_J": P_J,
        "P_ep": P_ep,
        "P_diff": P_diff,
        "P_esc": P_esc,
        "Q_ret": Q_ret,
        "du_dt": du_dt,
        "residual": du_dt - Q_ret,
    }


def _append_temporal_metrics(
    target: dict[str, list[np.ndarray]],
    *,
    interval: Mapping[str, np.ndarray],
    derived: Mapping[str, np.ndarray],
    center_indices: np.ndarray,
    volume_weights: np.ndarray,
    volume_total: float,
    xi_m: float,
    delta0_J: float,
) -> None:
    del center_indices
    target["time_ps"].append(np.asarray(interval["time_ps"], dtype=float))
    target["dt_ps"].append(np.asarray(interval["dt_ps"], dtype=float))
    for source, destination in (
        ("P_spec", "p99_abs_P_spec_W_m3"),
        ("P_delta", "p99_abs_P_delta_W_m3"),
        ("P_delta_cond", "p99_abs_P_delta_cond_W_m3"),
        ("P_q", "p99_abs_P_q_W_m3"),
        ("P_path", "p99_abs_P_path_W_m3"),
        ("Q_ret", "p99_abs_Q_ret_W_m3"),
        ("residual", "p99_abs_residual_W_m3"),
    ):
        target[destination].append(np.nanpercentile(np.abs(interval[source]), 99.0, axis=1))
    for source, destination, sign in (
        ("P_spec", "integrated_P_spec_W", 1.0),
        ("P_delta", "integrated_P_delta_W", 1.0),
        ("P_delta_cond", "integrated_P_delta_cond_W", 1.0),
        ("P_q", "integrated_P_q_W", 1.0),
        ("P_path", "integrated_P_path_W", 1.0),
        ("P_J", "integrated_P_J_W", 1.0),
        ("P_ep", "integrated_minus_P_ep_W", -1.0),
        ("P_diff", "integrated_P_diff_W", 1.0),
        ("Q_ret", "integrated_Q_ret_W", 1.0),
        ("residual", "integrated_residual_W", 1.0),
        ("P_esc", "integrated_P_esc_W", 1.0),
    ):
        target[destination].append(sign * np.sum(interval[source] * volume_weights[None, :], axis=1))
    u_e = np.asarray(derived["u_e_center"], dtype=float)
    u_cond = np.asarray(derived["u_cond_center"], dtype=float)
    u_ph = np.asarray(derived["u_ph_center"], dtype=float)
    for values, destination in (
        (u_e, "mean_u_e_J_m3"),
        (u_cond, "mean_u_cond_J_m3"),
        (u_e - u_cond, "mean_u_qp_J_m3"),
        (u_ph, "mean_u_ph_J_m3"),
    ):
        midpoint = 0.5 * (values[:-1] + values[1:])
        target[destination].append(
            np.sum(midpoint * volume_weights[None, :], axis=1) / max(volume_total, 1.0e-300)
        )
    delta_center = np.asarray(derived["delta_center"], dtype=float)
    q_center = np.asarray(derived["q_center"], dtype=float)
    Te_center = np.asarray(derived["Te_center"], dtype=float)
    Tph_center = np.asarray(derived["Tph_center"], dtype=float)
    target["min_delta_over_delta0"].append(np.nanmin(delta_center[1:], axis=1) / delta0_J)
    target["max_Te_K"].append(np.nanmax(Te_center[1:], axis=1))
    target["max_Tph_K"].append(np.nanmax(Tph_center[1:], axis=1))
    target["max_q_xi"].append(np.nanmax(q_center[1:], axis=1) * xi_m)
    clip = np.asarray(derived["clip_center"], dtype=bool)
    target["catalog_clipped_fraction"].append(np.mean(clip[1:], axis=1))
    for component in ("Te", "Tph", "delta", "q"):
        component_clip = np.asarray(derived[f"clip_{component}_center"], dtype=bool)
        target[f"catalog_{component}_clipped_fraction"].append(
            np.mean(component_clip[1:], axis=1)
        )
    strict_clip = np.asarray(derived["strict_q_clip_center"], dtype=bool)[1:]
    js_abs = np.asarray(derived["js_abs_center_A_m2"], dtype=float)[1:]
    js_reference = float(
        np.asarray(derived["js_reference_A_m2"], dtype=float).reshape(-1)[0]
    )
    target["strict_q_clipped_fraction"].append(np.mean(strict_clip, axis=1))
    conditional_median: list[float] = []
    conditional_p95: list[float] = []
    conditional_max: list[float] = []
    conditional_p95_ratio: list[float] = []
    conditional_near_zero: list[float] = []
    for row_mask, row_current in zip(strict_clip, js_abs):
        selected = np.abs(row_current[row_mask & np.isfinite(row_current)])
        if selected.size == 0:
            conditional_median.append(np.nan)
            conditional_p95.append(np.nan)
            conditional_max.append(np.nan)
            conditional_p95_ratio.append(np.nan)
            conditional_near_zero.append(np.nan)
            continue
        p95 = float(np.nanpercentile(selected, 95.0))
        conditional_median.append(float(np.nanmedian(selected)))
        conditional_p95.append(p95)
        conditional_max.append(float(np.nanmax(selected)))
        conditional_p95_ratio.append(p95 / max(js_reference, 1.0e-300))
        conditional_near_zero.append(
            float(np.mean(selected <= 1.0e-3 * max(js_reference, 1.0e-300)))
        )
    target["strict_q_clipped_js_median_A_m2"].append(
        np.asarray(conditional_median, dtype=float)
    )
    target["strict_q_clipped_js_p95_A_m2"].append(
        np.asarray(conditional_p95, dtype=float)
    )
    target["strict_q_clipped_js_max_A_m2"].append(
        np.asarray(conditional_max, dtype=float)
    )
    target["strict_q_clipped_js_p95_over_catalog_max"].append(
        np.asarray(conditional_p95_ratio, dtype=float)
    )
    target["strict_q_clipped_js_near_zero_fraction"].append(
        np.asarray(conditional_near_zero, dtype=float)
    )


def _full_state_maps(
    *,
    position: int,
    raw: Mapping[str, np.ndarray],
    derived: Mapping[str, np.ndarray],
    catalog: EnergyProjectionCatalog,
    xi_m: float,
) -> dict[str, np.ndarray]:
    Te = np.asarray(raw["Te_snapshot_K"], dtype=float)[position]
    Tph = np.asarray(raw["Tph_snapshot_K"], dtype=float)[position]
    delta = np.asarray(derived["delta_full"], dtype=float)[position]
    q_abs = np.asarray(derived["q_full"], dtype=float)[position]
    lookup = catalog.evaluate(Te_K=Te, Tph_K=Tph, delta_J=delta, q_abs_m_inv=q_abs)
    u_cond = catalog.condensation_energy(delta)
    return {
        "delta_over_delta0": delta / catalog.delta0_J,
        "q_xi": q_abs * float(xi_m),
        "Te_K": Te,
        "Tph_K": Tph,
        "u_e_J_m3": lookup["u_e_J_m3"],
        "u_cond_J_m3": u_cond,
        "u_qp_J_m3": lookup["u_e_J_m3"] - u_cond,
        "u_ph_J_m3": lookup["u_ph_J_m3"],
    }


def _selected_interval_maps(
    *,
    local_end: int,
    raw: Mapping[str, np.ndarray],
    derived: Mapping[str, np.ndarray],
    interval: Mapping[str, np.ndarray],
    catalog: EnergyProjectionCatalog,
    projector: _SnapshotProjector,
    central_mask: np.ndarray,
    sigma_n_S_m: float,
    bath_K: float,
    xi_m: float,
) -> dict[str, np.ndarray]:
    start = local_end - 1
    Te0 = np.asarray(raw["Te_snapshot_K"], dtype=float)[start]
    Te1 = np.asarray(raw["Te_snapshot_K"], dtype=float)[local_end]
    Tph0 = np.asarray(raw["Tph_snapshot_K"], dtype=float)[start]
    Tph1 = np.asarray(raw["Tph_snapshot_K"], dtype=float)[local_end]
    delta0 = np.asarray(derived["delta_full"], dtype=float)[start]
    delta1 = np.asarray(derived["delta_full"], dtype=float)[local_end]
    q0 = np.asarray(derived["q_full"], dtype=float)[start]
    q1 = np.asarray(derived["q_full"], dtype=float)[local_end]
    lookup0 = catalog.evaluate(Te_K=Te0, Tph_K=Tph0, delta_J=delta0, q_abs_m_inv=q0)
    lookup1 = catalog.evaluate(Te_K=Te1, Tph_K=Tph1, delta_J=delta1, q_abs_m_inv=q1)
    kappa0 = np.zeros_like(Te0)
    kappa1 = np.zeros_like(Te1)
    kappa0[central_mask] = lookup0["kappa_W_m_K"][central_mask]
    kappa1[central_mask] = lookup1["kappa_W_m_K"][central_mask]
    P_diff0 = projector.diffusion_power_density(Te0[None, :], kappa0[None, :], active_mask=central_mask, bath_K=bath_K)[0]
    P_diff1 = projector.diffusion_power_density(Te1[None, :], kappa1[None, :], active_mask=central_mask, bath_K=bath_K)[0]
    phi = np.asarray(raw["phi_snapshot_V"], dtype=float)
    P_J_pair = projector.joule_power_density(phi[start : local_end + 1], sigma_n_S_m=sigma_n_S_m)
    dt_s = float(np.asarray(interval["dt_ps"])[start]) * 1.0e-12
    u0 = lookup0["u_e_J_m3"]
    u1 = lookup1["u_e_J_m3"]
    u_delta = catalog.electronic_energy(Te0, delta1, q0)
    u_new_spectrum = catalog.electronic_energy(Te0, delta1, q1)
    u_q_first = catalog.electronic_energy(Te0, delta0, q1)
    P_delta = (u_delta - u0) / dt_s
    P_q = (u_new_spectrum - u_delta) / dt_s
    P_spec = (u_new_spectrum - u0) / dt_s
    P_path = P_delta - (u_new_spectrum - u_q_first) / dt_s
    u_cond0 = catalog.condensation_energy(delta0)
    u_cond1 = catalog.condensation_energy(delta1)
    P_delta_cond = (u_cond1 - u_cond0) / dt_s
    P_thermal = (u1 - u_new_spectrum) / dt_s
    P_CdT = 0.5 * (lookup0["C_e_J_m3_K"] + lookup1["C_e_J_m3_K"]) * (Te1 - Te0) / dt_s
    P_J = np.mean(P_J_pair, axis=0)
    P_ep = 0.5 * (lookup0["P_ep_W_m3"] + lookup1["P_ep_W_m3"])
    P_diff = 0.5 * (P_diff0 + P_diff1)
    Q_ret = P_J + P_diff - P_ep
    du_dt = (u1 - u0) / dt_s
    u_cond = u_cond1
    return {
        "delta_over_delta0": delta1 / catalog.delta0_J,
        "q_xi": q1 * float(xi_m),
        "Te_K": Te1,
        "Tph_K": Tph1,
        "u_e_J_m3": u1,
        "u_cond_J_m3": u_cond,
        "u_qp_J_m3": u1 - u_cond,
        "u_ph_J_m3": lookup1["u_ph_J_m3"],
        "P_delta_W_m3": P_delta,
        "P_delta_cond_W_m3": P_delta_cond,
        "P_delta_qp_W_m3": P_delta - P_delta_cond,
        "P_q_W_m3": P_q,
        "P_spec_W_m3": P_spec,
        "P_path_W_m3": P_path,
        "P_J_W_m3": P_J,
        "minus_P_ep_W_m3": -P_ep,
        "P_diff_W_m3": P_diff,
        "P_esc_W_m3": 0.5 * (lookup0["P_esc_W_m3"] + lookup1["P_esc_W_m3"]),
        "Q_ret_W_m3": Q_ret,
        "du_e_dt_W_m3": du_dt,
        "P_thermal_W_m3": P_thermal,
        "P_CdT_W_m3": P_CdT,
        "P_T_nonlinear_W_m3": P_thermal - P_CdT,
        "residual_W_m3": du_dt - Q_ret,
        "catalog_clipped": catalog.clipping_mask(Te1, Tph1, delta1, q1).astype(float),
    }


def _load_snapshot_times(path: str | Path) -> np.ndarray:
    with np.load(Path(path), allow_pickle=False) as data:
        if "snapshot_t_ps" in data.files:
            return np.asarray(data["snapshot_t_ps"], dtype=float).reshape(-1)
        return np.asarray(data["snapshot_t_s"], dtype=float).reshape(-1) / 1.0e-12


def _strictly_increasing_sample_indices(
    times_ps: np.ndarray,
    candidate_indices: np.ndarray,
    *,
    duplicate_tolerance_ps: float = 1.0e-9,
) -> np.ndarray:
    """Keep the later state when persisted snapshots repeat one timestamp."""

    times = np.asarray(times_ps, dtype=float).reshape(-1)
    candidates = np.asarray(candidate_indices, dtype=np.int64).reshape(-1)
    selected: list[int] = []
    for raw_index in candidates:
        index = int(raw_index)
        value = float(times[index])
        if not np.isfinite(value):
            raise ValueError(f"Snapshot time at index {index} is not finite.")
        if not selected:
            selected.append(index)
            continue
        previous_index = selected[-1]
        previous = float(times[previous_index])
        if value > previous:
            selected.append(index)
            continue
        if abs(value - previous) <= float(duplicate_tolerance_ps):
            # Persisted final payloads may repeat the requested terminal time.
            # The later row is the authoritative final state.
            selected[-1] = index
            continue
        raise ValueError(
            "Persisted snapshot times move backwards: "
            f"index {previous_index} has {previous:g} ps and index {index} "
            f"has {value:g} ps."
        )
    if len(selected) < 2:
        raise ValueError("Fewer than two strictly increasing snapshot times remain.")
    return np.asarray(selected, dtype=np.int64)


def _history_time_ps(history: Mapping[str, Any]) -> np.ndarray:
    if "t_ps" in history:
        return np.asarray(history["t_ps"], dtype=float).reshape(-1)
    return np.asarray(history.get("t_s", []), dtype=float).reshape(-1) / 1.0e-12


def _photon_time_ps(history: Mapping[str, Any]) -> float:
    time_ps = _history_time_ps(history)
    applied = np.asarray(history.get("photon_applied", []), dtype=bool).reshape(-1)
    count = min(time_ps.size, applied.size)
    hits = np.flatnonzero(applied[:count])
    return float(time_ps[hits[0]]) if hits.size else np.nan


def _read_exact(stream: Any, nbytes: int) -> bytes:
    blocks: list[bytes] = []
    remaining = int(nbytes)
    while remaining:
        block = stream.read(remaining)
        if not block:
            raise EOFError(f"Compressed NPY member ended {remaining} bytes early.")
        blocks.append(block)
        remaining -= len(block)
    return b"".join(blocks)


def _first_npz_array(data: Any, keys: Sequence[str]) -> np.ndarray:
    for key in keys:
        if key in data.files:
            return np.asarray(data[key], dtype=float)
    raise KeyError(f"None of the required NPZ arrays is present: {', '.join(keys)}")


__all__ = [
    "CompressedNpzRowStream",
    "EnergyProjectionCatalog",
    "SupercurrentProjectionCatalog",
    "extract_energy_projection_diagnostics",
]
