"""
Project configuration utilities for pySNSPD.

This module is intentionally limited to reading, normalizing and validating
configuration files. It does not create folders or write run data; that
responsibility belongs to ``pysnspd.io.manager``.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any, Mapping
import os


try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "pySNSPD requires PyYAML to read .yaml configuration files. "
        "Install it with: python -m pip install pyyaml"
    ) from exc


class ConfigError(ValueError):
    """Raised when a pySNSPD configuration file is missing required fields."""


def load_config(config_path: str | Path) -> dict[str, Any]:
    """
    Load a YAML configuration file.

    Parameters
    ----------
    config_path:
        Path to the YAML configuration file.

    Returns
    -------
    dict
        Configuration dictionary with a small ``_meta`` block containing the
        absolute path to the source file.
    """
    path = Path(config_path).expanduser().resolve()

    if not path.exists():
        raise ConfigError(f"Configuration file does not exist: {path}")

    if not path.is_file():
        raise ConfigError(f"Configuration path is not a file: {path}")

    with path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    if cfg is None:
        raise ConfigError(f"Configuration file is empty: {path}")

    if not isinstance(cfg, dict):
        raise ConfigError(f"Configuration file must contain a YAML mapping: {path}")

    cfg = deepcopy(cfg)
    cfg.setdefault("_meta", {})
    cfg["_meta"]["config_path"] = str(path)

    return cfg


def validate_config(
    config: Mapping[str, Any],
    *,
    require_big_data_root_exists: bool = True,
) -> dict[str, Any]:
    """
    Validate and normalize a pySNSPD configuration dictionary.

    The function mutates neither the input mapping nor the source YAML file.
    It returns a normalized copy. In particular, ``project.big_data_root`` is
    expanded to an absolute path.

    Parameters
    ----------
    config:
        Configuration dictionary, usually returned by :func:`load_config`.
    require_big_data_root_exists:
        If True, ``project.big_data_root`` must already exist.

    Returns
    -------
    dict
        Normalized configuration dictionary.

    Raises
    ------
    ConfigError
        If required fields are missing or invalid.
    """
    if not isinstance(config, Mapping):
        raise ConfigError("Configuration must be a mapping/dictionary.")

    cfg = deepcopy(dict(config))

    required_sections = [
        "project",
        "parallel",
        "material",
        "calibration",
        "bias",
        "mesh",
        "catalogs",
        "ss_run",
        "photon_run",
        "circuit",
    ]
    for section in required_sections:
        _require_section(cfg, section)

    project = _require_section(cfg, "project")
    _require_string(project, "name", "project")
    _require_string(project, "default_run_name", "project")

    big_data_root = _require_string(project, "big_data_root", "project")
    if "CHANGE/ME" in big_data_root or big_data_root.strip() == "":
        raise ConfigError(
            "project.big_data_root is not configured. "
            "Set it to an absolute external data folder, for example "
            "/home/jdiaz/scratch/big_data."
        )

    root = Path(os.path.expandvars(os.path.expanduser(big_data_root))).resolve()
    project["big_data_root"] = str(root)

    if require_big_data_root_exists:
        if not root.exists():
            raise ConfigError(f"project.big_data_root does not exist: {root}")
        if not root.is_dir():
            raise ConfigError(f"project.big_data_root is not a directory: {root}")
        if not os.access(root, os.W_OK):
            raise ConfigError(f"project.big_data_root is not writable: {root}")

    parallel = _require_section(cfg, "parallel")
    _require_bool(parallel, "enabled", "parallel")
    workers = _require_positive_int(parallel, "workers", "parallel")
    backend = _require_string(parallel, "backend", "parallel")
    if backend not in {"process", "thread", "serial"}:
        raise ConfigError(
            "parallel.backend must be one of: process, thread, serial. "
            f"Got: {backend}"
        )
    if backend == "serial" and workers != 1:
        raise ConfigError("parallel.backend='serial' requires parallel.workers=1.")

    material = _require_section(cfg, "material")
    _require_string(material, "name", "material")
    for key in [
        "Tc_K",
        "sigma_n_S_m",
        "lambda_L_m",
        "thickness_m",
        "width_m",
    ]:
        _require_positive_number(material, key, "material")

    if "D_m2_s" in material:
        _require_positive_number(material, "D_m2_s", "material")

    # gTDGL/KWT relaxation times at Tc.  These are optional for backward
    # compatibility, but when present they are validated here instead of being
    # silently accepted by the solver.  Values may be supplied either in seconds
    # or picoseconds; use exactly one alias per physical parameter.
    _validate_optional_time_parameter(
        material,
        section="material",
        label="tau_ee_Tc",
        names_s=("tau_ee_Tc_s", "tau_ee_s", "tau_ee_at_Tc_s"),
        names_ps=("tau_ee_Tc_ps", "tau_ee_ps", "tau_ee_at_Tc_ps"),
        default_s=5.0e-12,
    )
    _validate_optional_time_parameter(
        material,
        section="material",
        label="tau_ep_Tc",
        names_s=("tau_ep_Tc_s", "tau_ep_s", "tau_ep_at_Tc_s"),
        names_ps=("tau_ep_Tc_ps", "tau_ep_ps", "tau_ep_at_Tc_ps"),
        default_s=24.7e-12,
    )

    calibration = _require_section(cfg, "calibration")
    _require_positive_number(calibration, "Ic_target_A", "calibration")

    if "n_gamma_sweep" not in calibration:
        calibration["n_gamma_sweep"] = 160
    _require_positive_int(calibration, "n_gamma_sweep", "calibration")

    if "gamma_max_fraction" not in calibration:
        calibration["gamma_max_fraction"] = 0.80
    _require_positive_number(calibration, "gamma_max_fraction", "calibration")

    if "D_warn_min_m2_s" not in calibration:
        calibration["D_warn_min_m2_s"] = 5.0e-5
    _require_positive_number(calibration, "D_warn_min_m2_s", "calibration")

    if "D_warn_max_m2_s" not in calibration:
        calibration["D_warn_max_m2_s"] = 5.0e-4
    _require_positive_number(calibration, "D_warn_max_m2_s", "calibration")

    bias = _require_section(cfg, "bias")
    _require_positive_number(bias, "T_bias_K", "bias")
    _require_nonnegative_number(bias, "I_bias_A", "bias")

    mesh = _require_section(cfg, "mesh")
    mesh_type = _require_string(mesh, "type", "mesh")
    if mesh_type != "delaunay":
        raise ConfigError(f"Only mesh.type='delaunay' is supported for now. Got: {mesh_type}")
    _require_positive_number(mesh, "target_spacing_m", "mesh")
    _require_int(mesh, "seed", "mesh")

    if "length_m" in mesh:
        _require_positive_number(mesh, "length_m", "mesh")

    catalogs = _require_section(cfg, "catalogs")
    dos = _require_nested_section(catalogs, "dos", "catalogs")
    for key in ["n_delta", "n_q", "n_energy", "n_matsubara"]:
        _require_positive_int(dos, key, "catalogs.dos")

    phase_space = _require_nested_section(catalogs, "phase_space", "catalogs")
    for key in ["n_Te", "n_Tph", "n_delta", "n_q", "n_omega"]:
        _require_positive_int(phase_space, key, "catalogs.phase_space")

    ss_run = _require_section(cfg, "ss_run")
    # SS runs are controlled by physical time, not by a maximum step count.
    # ``max_steps`` is deliberately not required here.  Older YAML files are
    # tolerated, but the canonical key is ``total_time_ps``.
    if "total_time_ps" not in ss_run:
        if "physical_time_ps" in ss_run:
            ss_run["total_time_ps"] = _require_positive_number(ss_run, "physical_time_ps", "ss_run")
        elif "max_steps" in ss_run:
            # Backward compatibility only: convert old max_steps*dt_s to ps.
            # This keeps legacy configs loadable without keeping max_steps as
            # a control parameter in the SS pipeline.
            old_steps = _require_positive_int(ss_run, "max_steps", "ss_run")
            old_dt_s = _require_positive_number(ss_run, "dt_s", "ss_run")
            ss_run["total_time_ps"] = float(old_steps) * float(old_dt_s) * 1.0e12
        else:
            ss_run["total_time_ps"] = 20.0
    _require_positive_number(ss_run, "total_time_ps", "ss_run")
    _require_positive_number(ss_run, "dt_s", "ss_run")
    _require_positive_number(ss_run, "convergence_tol", "ss_run")
    if "snapshots" not in ss_run:
        ss_run["snapshots"] = 8
    _require_positive_int(ss_run, "snapshots", "ss_run")

    photon_run = _require_section(cfg, "photon_run")
    _require_positive_number(photon_run, "photon_wavelength_m", "photon_run")
    _require_positive_int(photon_run, "max_steps", "photon_run")
    _require_positive_number(photon_run, "dt_s", "photon_run")
    _require_positive_number(photon_run, "bubble_radius_m", "photon_run")

    circuit = _require_section(cfg, "circuit")
    _require_positive_number(circuit, "R_load_ohm", "circuit")
    _require_positive_number(circuit, "L_bias_H", "circuit")
    _require_positive_number(circuit, "C_rf_F", "circuit")

    return cfg


def summarize_config(config: Mapping[str, Any]) -> str:
    """
    Return a compact human-readable summary of the project configuration.
    """
    cfg = validate_config(config, require_big_data_root_exists=False)

    material_tau_ee_Tc_s = _validate_optional_time_parameter(
        cfg["material"],
        section="material",
        label="tau_ee_Tc",
        names_s=("tau_ee_Tc_s", "tau_ee_s", "tau_ee_at_Tc_s"),
        names_ps=("tau_ee_Tc_ps", "tau_ee_ps", "tau_ee_at_Tc_ps"),
        default_s=5.0e-12,
    )
    material_tau_ep_Tc_s = _validate_optional_time_parameter(
        cfg["material"],
        section="material",
        label="tau_ep_Tc",
        names_s=("tau_ep_Tc_s", "tau_ep_s", "tau_ep_at_Tc_s"),
        names_ps=("tau_ep_Tc_ps", "tau_ep_ps", "tau_ep_at_Tc_ps"),
        default_s=24.7e-12,
    )

    lines = [
        "pySNSPD project configuration",
        f"project.name             : {cfg['project']['name']}",
        f"project.big_data_root    : {cfg['project']['big_data_root']}",
        f"project.default_run_name : {cfg['project']['default_run_name']}",
        f"parallel.enabled         : {cfg['parallel']['enabled']}",
        f"parallel.workers         : {cfg['parallel']['workers']}",
        f"parallel.backend         : {cfg['parallel']['backend']}",
        f"material.name            : {cfg['material']['name']}",
        f"material.Tc_K            : {cfg['material']['Tc_K']}",
        f"material.width_m         : {cfg['material']['width_m']}",
        f"material.thickness_m     : {cfg['material']['thickness_m']}",
        f"material.tau_ee_Tc_ps    : {material_tau_ee_Tc_s / 1.0e-12}",
        f"material.tau_ep_Tc_ps    : {material_tau_ep_Tc_s / 1.0e-12}",
        f"bias.T_bias_K            : {cfg['bias']['T_bias_K']}",
        f"bias.I_bias_A            : {cfg['bias']['I_bias_A']}",
        f"calibration.Ic_target_A : {cfg['calibration']['Ic_target_A']}",
        f"calibration.jc_target_A_m2 : {cfg['calibration']['Ic_target_A'] / (cfg['material']['width_m'] * cfg['material']['thickness_m'])}",
        f"mesh.type                : {cfg['mesh']['type']}",
        f"mesh.target_spacing_m    : {cfg['mesh']['target_spacing_m']}",
        f"mesh.seed                : {cfg['mesh']['seed']}",
    ]
    return "\n".join(lines)


def _validate_optional_time_parameter(
    config: Mapping[str, Any],
    *,
    section: str,
    label: str,
    names_s: tuple[str, ...],
    names_ps: tuple[str, ...],
    default_s: float,
) -> float:
    """Validate one optional time parameter and return its value in seconds.

    The config supports both seconds and picoseconds aliases for convenience,
    but accepting multiple aliases for the same physical time is ambiguous.
    """
    found: list[tuple[str, float]] = []

    for name in names_s:
        if name in config:
            found.append((name, _require_positive_number(config, name, section)))

    for name in names_ps:
        if name in config:
            found.append((name, _require_positive_number(config, name, section) * 1.0e-12))

    if len(found) > 1:
        names = ", ".join(name for name, _ in found)
        raise ConfigError(
            f"{section}.{label} has multiple aliases configured ({names}). "
            f"Use exactly one of {names_s + names_ps}."
        )

    if found:
        return float(found[0][1])

    if default_s <= 0.0:
        raise ConfigError(f"Internal default for {section}.{label} must be positive.")
    return float(default_s)

def _require_section(config: Mapping[str, Any], section: str) -> dict[str, Any]:
    if section not in config:
        raise ConfigError(f"Missing required section: {section}")
    value = config[section]
    if not isinstance(value, dict):
        raise ConfigError(f"Section '{section}' must be a mapping/dictionary.")
    return value


def _require_nested_section(
    config: Mapping[str, Any],
    section: str,
    parent_name: str,
) -> dict[str, Any]:
    if section not in config:
        raise ConfigError(f"Missing required section: {parent_name}.{section}")
    value = config[section]
    if not isinstance(value, dict):
        raise ConfigError(f"Section '{parent_name}.{section}' must be a mapping/dictionary.")
    return value


def _require_string(config: Mapping[str, Any], key: str, section: str) -> str:
    if key not in config:
        raise ConfigError(f"Missing required key: {section}.{key}")
    value = config[key]
    if not isinstance(value, str):
        raise ConfigError(f"{section}.{key} must be a string.")
    if value.strip() == "":
        raise ConfigError(f"{section}.{key} must not be empty.")
    return value


def _require_bool(config: Mapping[str, Any], key: str, section: str) -> bool:
    if key not in config:
        raise ConfigError(f"Missing required key: {section}.{key}")
    value = config[key]
    if not isinstance(value, bool):
        raise ConfigError(f"{section}.{key} must be a boolean.")
    return value


def _require_positive_int(config: Mapping[str, Any], key: str, section: str) -> int:
    value = _require_int(config, key, section)
    if value <= 0:
        raise ConfigError(f"{section}.{key} must be positive.")
    return value


def _require_int(config: Mapping[str, Any], key: str, section: str) -> int:
    if key not in config:
        raise ConfigError(f"Missing required key: {section}.{key}")

    value = config[key]

    if isinstance(value, bool):
        raise ConfigError(f"{section}.{key} must be an integer.")

    if isinstance(value, int):
        return value

    if isinstance(value, str):
        try:
            parsed = int(value)
        except ValueError as exc:
            raise ConfigError(f"{section}.{key} must be an integer.") from exc

        if isinstance(config, dict):
            config[key] = parsed
        return parsed

    raise ConfigError(f"{section}.{key} must be an integer.")


def _require_positive_number(config: Mapping[str, Any], key: str, section: str) -> float:
    if key not in config:
        raise ConfigError(f"Missing required key: {section}.{key}")

    value = config[key]

    if isinstance(value, bool):
        raise ConfigError(f"{section}.{key} must be a number.")

    if isinstance(value, (int, float)):
        parsed = float(value)
    elif isinstance(value, str):
        try:
            parsed = float(value)
        except ValueError as exc:
            raise ConfigError(f"{section}.{key} must be a number.") from exc
    else:
        raise ConfigError(f"{section}.{key} must be a number.")

    if parsed <= 0.0:
        raise ConfigError(f"{section}.{key} must be positive.")

    if isinstance(config, dict):
        config[key] = parsed

    return parsed


def _require_nonnegative_number(config: Mapping[str, Any], key: str, section: str) -> float:
    if key not in config:
        raise ConfigError(f"Missing required key: {section}.{key}")

    value = config[key]

    if isinstance(value, bool):
        raise ConfigError(f"{section}.{key} must be a number.")

    if isinstance(value, (int, float)):
        parsed = float(value)
    elif isinstance(value, str):
        try:
            parsed = float(value)
        except ValueError as exc:
            raise ConfigError(f"{section}.{key} must be a number.") from exc
    else:
        raise ConfigError(f"{section}.{key} must be a number.")

    if parsed < 0.0:
        raise ConfigError(f"{section}.{key} must be nonnegative.")

    if isinstance(config, dict):
        config[key] = parsed

    return parsed
