#!/usr/bin/env python3
"""Patch pySNSPD SS-run pipeline for snapshot diagnostics and time-limited SS runs."""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path.cwd()


def patch_pipeline() -> None:
    path = ROOT / "pipelines/02_ss_run_template.py"
    text = path.read_text(encoding="utf-8")

    import_old = "from pysnspd.gtdgl.state_io import save_relaxation_history_npz, save_stationary_state_npz"
    import_new = (
        "from pysnspd.gtdgl.state_io import save_relaxation_history_npz, save_stationary_state_npz\n"
        "from pysnspd.gtdgl.snapshot_diagnostics import (\n"
        "    save_ss_snapshot_bundle_npz,\n"
        "    write_ss_snapshot_power_diagnostics,\n"
        ")"
    )
    if "pysnspd.gtdgl.snapshot_diagnostics" not in text:
        if import_old not in text:
            raise SystemExit(f"Could not find state_io import in {path}")
        text = text.replace(import_old, import_new)

    # Remove the runtime use of --ss-steps.  The parser may keep a deprecated
    # compatibility flag in some branches, but the solver is now governed only by
    # a physical time limit.
    time_old = re.compile(
        r"if args\.ss_time_ps is not None:\s*\n"
        r"\s*total_time_ps = float\(args\.ss_time_ps\)\s*\n"
        r"\s*elif args\.ss_steps is not None:\s*\n"
        r"\s*total_time_ps = float\(args\.ss_steps\) \* float\(args\.ss_dt_fs\) \* 1\.0e-3\s*\n"
        r"\s*else:\s*\n"
        r"\s*total_time_ps = 20\.0",
        re.M,
    )
    time_new = (
        "ss_run_cfg = cfg.get(\"ss_run\", {}) if isinstance(cfg, dict) else {}\n"
        "    if args.ss_time_ps is not None:\n"
        "        total_time_ps = float(args.ss_time_ps)\n"
        "    else:\n"
        "        total_time_ps = float(\n"
        "            ss_run_cfg.get(\"total_time_ps\", ss_run_cfg.get(\"physical_time_ps\", 20.0))\n"
        "        )"
    )
    text2, n = time_old.subn(time_new, text, count=1)
    if n == 0 and "ss_run_cfg = cfg.get(\"ss_run\"" not in text:
        raise SystemExit("Could not patch total_time_ps block in 02_ss_run_template.py")
    text = text2

    text = text.replace(
        "steps=None if args.ss_steps is None else int(args.ss_steps),",
        "steps=None,",
    )

    history_line = "history_npz = save_relaxation_history_npz(result.history, raw_ss / \"relaxation_history.npz\")"
    inserted = (
        history_line
        + "\n"
        + "    snapshots_npz = save_ss_snapshot_bundle_npz(result.history, raw_ss / \"stationary_snapshots.npz\")\n"
        + "    snapshot_power_npz = None\n"
        + "    power_table_path = raw_pre / \"power_table_catalog.npz\"\n"
        + "    if power_table_path.exists():\n"
        + "        snapshot_power_npz = write_ss_snapshot_power_diagnostics(\n"
        + "            history=result.history,\n"
        + "            state=result.state,\n"
        + "            power_table_npz=power_table_path,\n"
        + "            output_path=raw_ss / \"snapshot_power_energy_diagnostics.npz\",\n"
        + "            sigma_n_S_m=float(material.sigma_n_S_m),\n"
        + "        )"
    )
    if "stationary_snapshots.npz" not in text:
        if history_line not in text:
            raise SystemExit(f"Could not find history save line in {path}")
        text = text.replace(history_line, inserted)

    out_line = '"relaxation_history_npz": str(history_npz),'
    out_new = (
        out_line
        + "\n"
        + '            "stationary_snapshots_npz": str(snapshots_npz),\n'
        + '            "snapshot_power_energy_diagnostics_npz": (\n'
        + '                str(snapshot_power_npz) if snapshot_power_npz is not None else None\n'
        + '            ),'
    )
    if "snapshot_power_energy_diagnostics_npz" not in text:
        if out_line not in text:
            raise SystemExit(f"Could not find outputs history line in {path}")
        text = text.replace(out_line, out_new, 1)

    print_line = 'print(f" relaxation_history_npz:{history_npz}")'
    print_new = (
        print_line
        + "\n"
        + '    print(f" stationary_snapshots_npz: {snapshots_npz}")\n'
        + '    if snapshot_power_npz is not None:\n'
        + '        print(f" snapshot_power_energy_diagnostics_npz: {snapshot_power_npz}")'
    )
    if "stationary_snapshots_npz:" not in text and print_line in text:
        text = text.replace(print_line, print_new, 1)

    path.write_text(text, encoding="utf-8")
    print(f"patched {path}")


def patch_config() -> None:
    path = ROOT / "configs/geminga_local_v3.yaml"
    if not path.exists():
        print(f"skip missing {path}")
        return
    text = path.read_text(encoding="utf-8")
    new_section = (
        "ss_run:\n"
        "  total_time_ps: 20.0\n"
        "  dt_s: 1.0e-15\n"
        "  convergence_tol: 1.0e-7\n"
        "  snapshots: 8\n"
    )
    pattern = re.compile(r"ss_run:\n(?:[ \t]+.*\n)+?(?=photon_run:)", re.M)
    text2, n = pattern.subn(new_section, text, count=1)
    if n == 0:
        print(f"warning: could not patch ss_run section in {path}; leaving it unchanged")
    else:
        path.write_text(text2, encoding="utf-8")
        print(f"patched {path}")


def main() -> int:
    patch_pipeline()
    patch_config()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
