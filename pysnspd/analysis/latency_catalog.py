"""Read-only, source-authenticated latency catalog for the Rs608 paper family.

Only accepted native 2D observations define tau_d. The separate display helper
can reflect missing positions without adding simulated samples to the catalog.
"""
from __future__ import annotations

from functools import lru_cache
import csv
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

from pysnspd.experimental.trigger_window import TriggerWindow, TriggerWindowConfig

EV_J = 1.602176634e-19
BASE_CAMPAIGN = "calibration_eta090_ee5_rs425_rs608_20261021_01"
BASE_NAMES = tuple(f"rs608_bias{i}_{w}" for i in ("15p5", "21p5") for w in (1550, 775))
CURRENTS_UA = (15.5, 18.5, 21.5)
# Changes in these fields do not change the native 2D constitutive experiment.
ALLOWED_VARIATIONS = frozenset(("current_uA", "wavelength_nm", "retained_fraction",
    "impact_y_nm", "deposited_energy_eV", "source_boundary", "threads", "physical_cores",
    "latency_only", "final_time_ns", "max_1d_wall_hours", "ceiling_ps",
    "no_trigger_time_ceiling_s"))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _equal(a, b):
    if isinstance(a, (float, int)) and isinstance(b, (float, int)):
        return bool(np.isclose(a, b, rtol=1e-12, atol=0))
    return a == b


def family_mismatches(manifest, reference, currents=CURRENTS_UA):
    """Match all inherited 2D inputs and every frozen baseline source pin."""
    if manifest.get("execution_mode", "normal") != "normal":
        return ["execution_mode:" + str(manifest["execution_mode"])]
    p = manifest.get("parameter_identity", {})
    mismatches = [f"parameter:{k}" for k, v in reference["parameter_identity"].items()
                  if k not in ALLOWED_VARIATIONS and not _equal(v, p.get(k))]
    mismatches += [f"unexpected parameter:{k}" for k in p
                   if k not in reference["parameter_identity"] and k not in ALLOWED_VARIATIONS]
    pins = manifest.get("source_pins_sha256", {})
    mismatches += [f"source:{k}" for k, v in reference["source_pins_sha256"].items()
                   if pins.get(k) != v]
    if currents is not None and p.get("current_uA") not in currents:
        mismatches.append("current outside requested set")
    elif currents is None and (not isinstance(p.get("current_uA"), (int, float))
                              or not np.isfinite(p["current_uA"]) or p["current_uA"] <= 0):
        mismatches.append("invalid positive finite current")
    return mismatches


def replay_observation(trace, saved):
    """Replay the existing confirmation/hysteresis controller, never a peak fit."""
    if len(trace) < 2:
        raise ValueError("Insufficient accepted observations")
    t = np.asarray(trace["time_ps"], float) * 1e-12
    v = np.asarray(trace["Vout_V"], float)
    if not np.isfinite(t).all() or not np.isfinite(v).all() or np.any(np.diff(t) <= 0):
        raise ValueError("Nonfinite or nonmonotone accepted observations")
    config = TriggerWindowConfig(**saved["config"])
    if abs(t[0] - config.event_time_s) > 1e-18:
        raise ValueError("Accepted observations must include the full event-time origin")
    controller = TriggerWindow(config)
    for time_s, voltage in zip(t, v):
        controller.observe_accepted(float(time_s), float(voltage))
        if controller.terminal:
            break
    replay = controller.summary()
    for key in ("status", "detected", "censored", "margin_complete"):
        if replay[key] != saved[key]:
            raise ValueError(f"Trigger replay mismatch: {key}")
    for key in ("crossing_time_s", "latency_s", "confirmation_time_s", "actual_stop_s"):
        a, b = replay[key], saved[key]
        if (a is None) != (b is None) or (a is not None and abs(a - b) > 1e-18):
            raise ValueError(f"Trigger replay mismatch: {key}")
    if saved["detected"] and not saved["margin_complete"]:
        raise ValueError("Incomplete post-trigger observation margin")
    return replay


def deduplicate(rows):
    """Reconcile censored extensions; conflicting physical observations fail closed."""
    unique = {}
    aliases = []
    for row in sorted(rows, key=lambda r: (not r["is_baseline"], r["run_dir"])):
        key = (row["current_uA"], round(row["retained_energy_eV"], 10), round(row["y0_nm"], 8))
        if key not in unique:
            unique[key] = row
            continue
        old = unique[key]
        if old["censored"] or row["censored"]:
            if not old.get("initial_state_sha256") or old.get("initial_state_sha256") != row.get("initial_state_sha256"):
                raise ValueError(f"Conflicting duplicate physical input (unmatched preparation): {key}")
            if old["censored"] and row["censored"]:
                winner = max((old, row), key=lambda r: r["censored_lower_bound_ps"])
                reason = "longest authenticated censored observation"
            else:
                winner, limit = (row, old) if old["censored"] else (old, row)
                if winner["tau_d_ps"] < limit["censored_lower_bound_ps"] - 1e-6:
                    raise ValueError(f"Conflicting duplicate physical input (crossing before censor bound): {key}")
                reason = "authenticated detection resolves an earlier censored observation"
            loser = row if winner is old else old
            # Redirect earlier alias receipts when a later extension wins.
            for alias in aliases:
                if alias["canonical"] == loser["run_dir"]:
                    alias["canonical"] = winner["run_dir"]
            unique[key] = winner
            aliases.append(dict(run_dir=loser["run_dir"], canonical=winner["run_dir"], reason=reason))
            continue
        if abs(old["tau_d_ps"] - row["tau_d_ps"]) > 1e-6:
            raise ValueError(f"Conflicting duplicate physical input: {key}")
        aliases.append({"run_dir": row["run_dir"], "canonical": old["run_dir"], "reason": "duplicate physical input"})
    return sorted(unique.values(), key=lambda r: (r["current_uA"], r["retained_energy_eV"], r["y0_nm"])), aliases


def load_verified_catalog(csv_path):
    """Reuse a frozen publication table, checking its receipt and typed rows.

    This reads a previously authenticated data product; it does not certify the
    current bytes of its native archives. Current manifest/source identity is
    checked separately before each reuse, and changed/new cases use full replay.
    """
    csv_path = Path(csv_path).resolve()
    audit_path = csv_path.with_name("audit.json")
    audit = json.loads(audit_path.read_text())
    if audit.get("schema") != "pysnspd.paper_latency.v1":
        raise ValueError("Unrecognized prior latency catalog receipt")
    if audit.get("output_sha256", {}).get(csv_path.name) != sha256(csv_path):
        raise ValueError("Prior latency catalog digest mismatch")
    bool_fields = {"is_baseline", "censored", "trigger_replayed"}
    int_fields = {"accepted_trace_rows"}
    float_fields = {"current_uA", "retained_energy_eV", "y0_nm", "tau_d_ps",
                    "censored_lower_bound_ps", "threshold_mV", "confirmation_ps",
                    "max_global_energy_fraction", "observation_horizon_ps"}
    rows = {}
    with csv_path.open(newline="") as stream:
        for record in csv.DictReader(stream):
            row = {}
            for key, value in record.items():
                if value == "":
                    row[key] = None
                elif key in bool_fields:
                    if value not in ("True", "False"):
                        raise ValueError("Invalid prior catalog boolean")
                    row[key] = value == "True"
                elif key in int_fields:
                    row[key] = int(value)
                elif key in float_fields:
                    row[key] = float(value)
                    if not np.isfinite(row[key]):
                        raise ValueError("Nonfinite prior catalog quantity")
                else:
                    row[key] = value
            if not row.get("trigger_replayed") or not row.get("manifest_sha256"):
                raise ValueError("Prior row has no authenticated trigger receipt")
            if row["run_dir"] in rows:
                raise ValueError("Duplicate prior catalog path")
            rows[row["run_dir"]] = row
    if len(rows) != audit.get("unique_configurations"):
        raise ValueError("Prior catalog row count differs from receipt")
    return rows, audit, dict(catalog=str(csv_path), catalog_sha256=sha256(csv_path),
        audit=str(audit_path), audit_sha256=sha256(audit_path),
        prior_generated_utc=audit.get("generated_utc"), reused_runs=[])


def observation_identity(config):
    """An observation horizon is censoring, not a change to the trigger physics."""
    return {k: v for k, v in config.items() if k != "no_trigger_time_ceiling_s"}


def catalog(database, source_root, reuse_catalog=None, currents=CURRENTS_UA):
    database, source_root = Path(database).resolve(), Path(source_root).resolve()
    anchors = [database / BASE_CAMPAIGN / "runs" / name for name in BASE_NAMES]
    reference = json.loads((anchors[0] / "manifest.json").read_text())
    for path in anchors:
        if family_mismatches(json.loads((path / "manifest.json").read_text()), reference):
            raise ValueError(f"Baseline family is internally inconsistent: {path}")

    prior_rows, prior_receipt = {}, None
    if reuse_catalog is not None:
        prior_rows, prior_audit, prior_receipt = load_verified_catalog(reuse_catalog)
        if (prior_audit.get("baseline_source_pins") != reference["source_pins_sha256"]
                or observation_identity(prior_audit.get("trigger_config", {})) != observation_identity(reference["observation"]["config"])):
            raise ValueError("Prior catalog is not the current four-baseline experiment")

    @lru_cache(maxsize=None)
    def source_digest(relative):
        path = (source_root / relative).resolve()
        if not path.is_relative_to(source_root):
            raise ValueError("Source pin leaves repository")
        return sha256(path)

    reference_config = reference["observation"]["config"]
    rows, rejected, seen_paths = [], [], set()
    # Recursive discovery includes future campaign layouts; realpath removes aliases.
    for path in sorted(database.rglob("manifest.json")):
        run = path.parent.resolve()
        if run in seen_paths:
            continue
        seen_paths.add(run)
        try:
            raw = path.read_bytes()
            m = json.loads(raw)
            if "parameter_identity" not in m or "source_pins_sha256" not in m:
                continue
            mismatch = family_mismatches(m, reference, currents=currents)
            if mismatch:
                rejected.append({"run_dir": str(run), "reason": "different family", "details": mismatch})
                continue
            extension = m.get("latency_observation_extension")
            if m.get("phase") == "0d" or (m.get("continuation") and not extension):
                rejected.append({"run_dir": str(run), "reason": "linked continuation, not independent input"})
                continue
            if extension:
                # Only a full-observation 2D extension can resolve a censored
                # event. Thermal/electrical continuations are not latency data.
                parent = Path(extension["parent_run_dir"]).resolve()
                parent_manifest = json.loads((parent / "manifest.json").read_text())
                if (m.get("phase") != "2d" or
                        extension.get("parent_manifest_sha256") != sha256(parent / "manifest.json") or
                        sha256(run / "source_initial_state.npz") != sha256(parent / "source_initial_state.npz") or
                        family_mismatches(parent_manifest, reference, currents=currents)):
                    raise ValueError("Unverified full-observation latency extension")
            if m.get("status") in ("RUNNING", "PREPARATION"):
                rejected.append({"run_dir": str(run), "reason": "live/incomplete", "status": m.get("status")})
                continue
            obs = m.get("observation", {})
            if obs.get("status") not in ("TRIGGER_MARGIN_COMPLETE", "TIME_CENSORED_NO_TRIGGER"):
                rejected.append({"run_dir": str(run), "reason": "no admitted terminal trigger window", "status": m.get("status")})
                continue
            if observation_identity(obs["config"]) != observation_identity(reference_config) or obs.get("observable") != "V_out_V":
                raise ValueError("Trigger definition differs from four-case reference")
            for relative, expected in m["source_pins_sha256"].items():
                if source_digest(relative) != expected:
                    raise ValueError(f"Current source differs from saved producer: {relative}")
            previous = prior_rows.get(str(run))
            if previous is not None and hashlib.sha256(raw).hexdigest() == previous["manifest_sha256"]:
                # The immutable table is the plotted input; no native archive read
                # is represented as a fresh replay or a fresh archive integrity check.
                rows.append(dict(previous, observation_horizon_ps=obs["config"]["no_trigger_time_ceiling_s"] * 1e12))
                prior_receipt["reused_runs"].append(str(run))
                continue
            if sha256(run / "checkpoint_latest.npz") != m["checkpoint_sha256"]:
                raise ValueError("Checkpoint digest mismatch")
            s, p = m["source"], m["parameter_identity"]
            energy = float(s["retained_energy_J"])
            if not 0 < energy <= float(s["photon_energy_J"]) * (1 + 1e-12):
                raise ValueError("Invalid retained photon energy")
            if not np.isclose(energy, s["photon_energy_J"] * p["retained_fraction"], rtol=1e-10, atol=0):
                raise ValueError("Retained energy and input fraction disagree")
            y0 = float(s.get("impact_y_nm", p.get("impact_y_nm", 0)))
            if not np.isclose(y0, p.get("impact_y_nm", 0), atol=1e-8) or abs(y0) > 40 + 1e-8:
                raise ValueError("Invalid transverse position")
            if abs(s.get("relative_total_energy_defect", 0)) > 1e-10:
                raise ValueError("Initial deposition fails its energy normalization")
            count = int(m["archive_counts"]["2d"]["trace"])
            with h5py.File(run / "data.h5", "r", swmr=True) as h:
                q = h["2d"]
                if len(q["trace"]) < count:
                    raise ValueError("Archive shorter than committed trace")
                trace = q["trace"][:count]
                if extension:
                    parent_count = int(parent_manifest["archive_counts"]["2d"]["trace"])
                    with h5py.File(parent / "data.h5", "r", swmr=True) as parent_archive:
                        parent_trace = parent_archive["2d/trace"][:parent_count]
                    if count < parent_count or trace[:parent_count].tobytes() != parent_trace.tobytes():
                        raise ValueError("Latency extension does not preserve the complete parent trace prefix")
                if not np.isclose(q.attrs["R_load_ohm"], 50):
                    raise ValueError("Readout load differs from reference")
                if not np.isclose(q["coordinates_m"][:, 1].max() - q["coordinates_m"][:, 1].min(), 80e-9, rtol=1e-10, atol=1e-15):
                    raise ValueError("Strip width differs from reference")
            replay = replay_observation(trace, obs)
            if path.read_bytes() != raw:
                raise ValueError("Manifest changed while extracting")
            rows.append(dict(case=run.name, run_dir=str(run), model_version=m["model_version"],
                is_baseline=run in anchors, current_uA=float(p["current_uA"]),
                retained_energy_eV=energy / EV_J, y0_nm=y0,
                tau_d_ps=None if replay["latency_s"] is None else replay["latency_s"] * 1e12,
                censored=replay["censored"], censored_lower_bound_ps=None if not replay["censored"] else
                    replay["censored_lower_bound_s"] * 1e12,
                status=m["status"], observation_status=obs["status"], accepted_trace_rows=count,
                threshold_mV=obs["config"]["threshold_V"] * 1e3,
                confirmation_ps=obs["config"]["confirmation_s"] * 1e12,
                observation_horizon_ps=obs["config"]["no_trigger_time_ceiling_s"] * 1e12,
                manifest_sha256=hashlib.sha256(raw).hexdigest(),
                checkpoint_sha256=m["checkpoint_sha256"],
                initial_state_sha256=sha256(run / "source_initial_state.npz"),
                accepted_trace_sha256=hashlib.sha256(trace.tobytes()).hexdigest(),
                max_global_energy_fraction=m.get("max_global_energy_fraction"),
                family_anchor=str(anchors[0]), trigger_replayed=True,
                physical_scope="Effective-model deterministic latency; not experimental validation"))
        except (ValueError, KeyError, OSError, TypeError) as exc:
            rejected.append({"run_dir": str(run), "reason": "failed authentication", "details": str(exc)})
    unique, aliases = deduplicate(rows)
    present = {r["run_dir"] for r in rows}
    missing = [str(p) for p in anchors if str(p) not in present]
    if missing:
        raise ValueError(f"Required original baselines did not authenticate: {missing}")
    return unique, dict(rejected=rejected, aliases=aliases, baseline_runs=[str(p) for p in anchors],
        baseline_source_pins=reference["source_pins_sha256"],
        trigger_config=reference_config, discovered_real_runs=len(seen_paths),
        prior_verified_catalog=prior_receipt,
        verification_scope=("Reused immutable previously authenticated table rows with unchanged current manifests and source pins; new/changed cases fully replayed"
                            if prior_receipt else "All rows authenticated from native archives in this extraction"))


def select_slices(rows, position_nm=0., energy_eV=None):
    """Use one declared position/energy for all three currents, never marginalize."""
    energy_rows = [r for r in rows if np.isclose(r["y0_nm"], position_nm, rtol=0, atol=1e-8)]
    if energy_eV is None:
        energies = sorted({round(r["retained_energy_eV"], 10) for r in rows})
        # Cover detected currents first; censored-only coverage cannot form a curve.
        def score(e):
            group = [r for r in rows if np.isclose(r["retained_energy_eV"], e, rtol=0, atol=1e-9)]
            return len({r["current_uA"] for r in group if not r["censored"]}), len(group), -e
        best = max(energies, key=score)
        energy_eV = next(r["retained_energy_eV"] for r in rows if abs(r["retained_energy_eV"] - best) < 1e-9)
    position_rows = [r for r in rows if np.isclose(r["retained_energy_eV"], energy_eV, rtol=0, atol=1e-9)]
    for name, group in (("energy", energy_rows), ("position", position_rows)):
        missing = set(CURRENTS_UA) - {r["current_uA"] for r in group if not r["censored"]}
        if missing:
            raise ValueError(f"{name} slice has no detected point for currents {sorted(missing)}")
    return energy_rows, position_rows, float(energy_eV)


def position_display_rows(rows, reflect=True):
    """Reflect missing positions for display only; retain measured partners."""
    display = [dict(r, symmetry_mirrored=False, source_y0_nm=r["y0_nm"]) for r in rows]
    if reflect:
        for r in rows:
            y = r["y0_nm"]
            if abs(y) < 1e-8:
                continue
            partner = any(s["current_uA"] == r["current_uA"] and
                np.isclose(s["retained_energy_eV"], r["retained_energy_eV"], rtol=0, atol=1e-9) and
                np.isclose(s["y0_nm"], -y, rtol=0, atol=1e-8) for s in rows)
            if not partner:
                display.append(dict(r, y0_nm=-y, symmetry_mirrored=True, source_y0_nm=y))
    return sorted(display, key=lambda r: (r["current_uA"], r["y0_nm"]))
