"""Coverage, stationarity, and regime summary for Z2 current sweeps."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
import yaml

from pysnspd.plotting.style import THESIS_DPI, THESIS_WIDTH_IN, apply_thesis_style


apply_thesis_style()

_STATUS_COLORS = {True: "tab:green", False: "tab:red", None: "0.72"}


def build_current_sweep_regime_summary(
    points: Sequence[Mapping[str, Any]],
    skipped: Sequence[Mapping[str, Any]],
    *,
    ohmic_relative_tolerance: float,
) -> dict[str, Any]:
    """Build a serializable summary without interpolating across missing runs."""
    cases = _merged_cases(points, skipped)
    fields = {
        "complete": "complete",
        "strict_stationary": "strict_stationarity_passes",
        "dynamic_stationary": "dynamic_stationarity_passes",
        "photon_ready": "photon_ready",
        "approximately_ohmic": "approximately_ohmic",
    }
    ranges = {
        name: _sampled_range(cases, field)
        for name, field in fields.items()
    }
    incomplete = [case for case in cases if case.get("complete") is not True]
    return {
        "schema_version": 1,
        "criteria": {
            "photon_ready": (
                "Stored final photon-ready gate, or the same final gate reconstructed "
                "from strict-or-dynamic mesoscopic stationarity plus contact, continuity, "
                "thermal, circuit, and phase-drive gates."
            ),
            "approximately_ohmic": (
                "Both |V_TDGL|/V_N,TDGL and |V_terminal|/V_N,terminal lie within "
                f"{100.0 * float(ohmic_relative_tolerance):.3g}% of unity."
            ),
            "range_policy": (
                "Ranges describe only sampled currents; no state is inferred between "
                "missing or incomplete cases."
            ),
        },
        "counts": {
            "discovered": int(len(cases)),
            "complete": int(sum(case.get("complete") is True for case in cases)),
            "incomplete": int(len(incomplete)),
            "strict_stationary": int(sum(case.get("strict_stationarity_passes") is True for case in cases)),
            "dynamic_stationary": int(sum(case.get("dynamic_stationarity_passes") is True for case in cases)),
            "photon_ready": int(sum(case.get("photon_ready") is True for case in cases)),
            "approximately_ohmic": int(sum(case.get("approximately_ohmic") is True for case in cases)),
        },
        "sampled_ranges": ranges,
        "incomplete_currents_uA": [
            float(case["current_uA"])
            for case in incomplete
            if np.isfinite(float(case.get("current_uA", np.nan)))
        ],
        "cases": [_serializable_case(case) for case in cases],
    }


def write_current_sweep_regime_summary(
    summary: Mapping[str, Any],
    output_path: str | Path,
) -> Path:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(dict(summary), stream, sort_keys=False, allow_unicode=True)
    return path


def plot_current_sweep_regime_summary(
    points: Sequence[Mapping[str, Any]],
    skipped: Sequence[Mapping[str, Any]],
    output_path: str | Path,
    *,
    ohmic_relative_tolerance: float,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot four compact panels for coverage, regimes, and solver cost."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    cases = _merged_cases(points, skipped)
    completed = [case for case in cases if case.get("complete") is True]
    currents = np.asarray(
        [float(case.get("current_uA", np.nan)) for case in cases], dtype=float
    )
    finite_currents = currents[np.isfinite(currents)]

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(THESIS_WIDTH_IN, 1.25 * THESIS_WIDTH_IN),
        sharex="col",
        constrained_layout=False,
    )
    fig.subplots_adjust(
        left=0.105,
        right=0.91,
        bottom=0.155,
        top=0.965,
        hspace=0.34,
        wspace=0.38,
    )

    _plot_status_matrix(axes[0, 0], cases)
    _plot_voltage_ratios(
        axes[0, 1],
        completed,
        ohmic_relative_tolerance=float(ohmic_relative_tolerance),
    )
    _plot_physical_regime(axes[1, 0], completed)
    _plot_integrator_cost(axes[1, 1], completed, cases)

    if finite_currents.size:
        pad = max(0.04 * float(np.ptp(finite_currents)), 1.5)
        limits = (float(np.min(finite_currents) - pad), float(np.max(finite_currents) + pad))
        for ax in axes.reshape(-1):
            ax.set_xlim(*limits)
    axes[1, 0].set_xlabel(r"$I_{\mathrm{bias}}$ [$\mu$A]")
    axes[1, 1].set_xlabel(r"$I_{\mathrm{bias}}$ [$\mu$A]")

    summary = build_current_sweep_regime_summary(
        points,
        skipped,
        ohmic_relative_tolerance=float(ohmic_relative_tolerance),
    )
    fig.text(
        0.5,
        0.072,
        _summary_footer(summary),
        ha="center",
        va="center",
        fontsize=7.2,
        linespacing=1.25,
    )
    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    return output


def _plot_status_matrix(ax: plt.Axes, cases: Sequence[Mapping[str, Any]]) -> None:
    rows = (
        ("complete", "Complete"),
        ("strict_stationarity_passes", "Strict SS"),
        ("dynamic_stationarity_passes", "Dynamic SS"),
        ("photon_ready", "Photon ready"),
    )
    for row, (field, _) in enumerate(rows):
        for case in cases:
            current = float(case.get("current_uA", np.nan))
            if not np.isfinite(current):
                continue
            status = _optional_bool(case.get(field))
            ax.scatter(
                current,
                row,
                marker="s",
                s=42.0,
                color=_STATUS_COLORS[status],
                edgecolor="0.2",
                linewidth=0.35,
                zorder=3,
            )
    ax.set_yticks(np.arange(len(rows)), [label for _, label in rows])
    ax.set_ylim(len(rows) - 0.45, -0.55)
    ax.set_title("Coverage and final gates")
    ax.grid(True, axis="x", linewidth=0.4, alpha=0.3)
    handles = [
        plt.Line2D([], [], marker="s", linestyle="", color=color, markeredgecolor="0.2", markersize=5.2, label=label)
        for color, label in (("tab:green", "Pass"), ("tab:red", "Fail"), ("0.72", "Unavailable"))
    ]
    ax.legend(handles=handles, loc="lower left", frameon=True, ncol=3, fontsize=6.7)


def _plot_voltage_ratios(
    ax: plt.Axes,
    completed: Sequence[Mapping[str, Any]],
    *,
    ohmic_relative_tolerance: float,
) -> None:
    x = _case_array(completed, "current_uA")
    central = _case_array(completed, "normal_voltage_ratio")
    terminal = _case_array(completed, "terminal_normal_voltage_ratio")
    tol = float(ohmic_relative_tolerance)
    ax.axhspan(1.0 - tol, 1.0 + tol, color="tab:green", alpha=0.12, label=rf"Ohmic $\pm${100*tol:.0f}%")
    ax.axhline(1.0, color="0.3", linewidth=0.8)
    ax.plot(x, central, "o-", color="tab:blue", label=r"$|V_{\mathrm{TDGL}}|/V_{N}$")
    ax.plot(x, terminal, "s--", color="tab:orange", label=r"$|V_{\mathrm{terminal}}|/V_{N}$")
    ax.set_ylabel("Normalized voltage")
    ax.set_title("Electrical approach to normal state")
    ax.grid(True, linewidth=0.4, alpha=0.3)
    ax.legend(loc="best", frameon=True, fontsize=6.8)


def _plot_physical_regime(ax: plt.Axes, completed: Sequence[Mapping[str, Any]]) -> None:
    x = _case_array(completed, "current_uA")
    mean_delta = _case_array(completed, "mean_delta_over_delta0")
    normal_like = _case_array(completed, "normal_like_fraction_final")
    ax.plot(x, mean_delta, "o-", color="tab:purple", label=r"Mean $|\Delta|/\Delta_0$")
    ax.plot(x, normal_like, "s--", color="tab:brown", label="Normal-like node fraction")
    ax.set_ylabel("Fraction / normalized gap")
    ax.set_ylim(0.0, 1.02)
    ax.set_title("Final condensate morphology")
    ax.grid(True, linewidth=0.4, alpha=0.3)
    ax.legend(loc="best", frameon=True, fontsize=6.8)


def _plot_integrator_cost(
    ax: plt.Axes,
    completed: Sequence[Mapping[str, Any]],
    cases: Sequence[Mapping[str, Any]],
) -> None:
    x = _case_array(completed, "current_uA")
    accepted = _case_array(completed, "accepted_steps")
    rejected_ratio = _case_array(completed, "rejected_over_accepted")
    ax.plot(x, accepted, "o-", color="tab:cyan", label="Accepted steps")
    ax.set_yscale("log")
    ax.set_ylabel("Accepted steps", color="tab:cyan")
    ax.tick_params(axis="y", colors="tab:cyan")
    twin = ax.twinx()
    twin.plot(x, rejected_ratio, "s--", color="tab:red", label="Rejected / accepted")
    twin.set_ylabel("Rejected / accepted", color="tab:red")
    twin.tick_params(axis="y", colors="tab:red")
    incomplete_x = [
        float(case.get("current_uA", np.nan))
        for case in cases
        if case.get("complete") is not True and np.isfinite(float(case.get("current_uA", np.nan)))
    ]
    if incomplete_x:
        ymin, _ = ax.get_ylim()
        ax.scatter(incomplete_x, np.full(len(incomplete_x), ymin), marker="x", color="0.35", s=24, label="Incomplete")
    ax.set_title("Integrator cost and failed coverage")
    ax.grid(True, linewidth=0.4, alpha=0.3)
    handles = ax.get_lines() + twin.get_lines()
    if incomplete_x:
        handles.append(ax.collections[-1])
    ax.legend(handles=handles, labels=[item.get_label() for item in handles], loc="best", frameon=True, fontsize=6.8)


def _merged_cases(
    points: Sequence[Mapping[str, Any]],
    skipped: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for point in points:
        if str(point.get("run_name", "")) == "synthetic_origin":
            continue
        case = dict(point)
        case.setdefault("complete", True)
        cases.append(case)
    for item in skipped:
        case = dict(item)
        case.setdefault("complete", False)
        cases.append(case)
    cases.sort(key=lambda item: (float(item.get("current_uA", np.inf)), str(item.get("run_name", ""))))
    return cases


def _sampled_range(cases: Sequence[Mapping[str, Any]], field: str) -> dict[str, Any]:
    currents = [
        float(case["current_uA"])
        for case in cases
        if case.get(field) is True and np.isfinite(float(case.get("current_uA", np.nan)))
    ]
    return {
        "currents_uA": currents,
        "minimum_uA": min(currents) if currents else None,
        "maximum_uA": max(currents) if currents else None,
        "n_sampled": int(len(currents)),
    }


def _summary_footer(summary: Mapping[str, Any]) -> str:
    ranges = summary.get("sampled_ranges", {})
    counts = summary.get("counts", {})
    photon = _compact_range(ranges.get("photon_ready", {}))
    dynamic = _compact_range(ranges.get("dynamic_stationary", {}))
    ohmic = _compact_range(ranges.get("approximately_ohmic", {}))
    return (
        f"Completed {counts.get('complete', 0)}/{counts.get('discovered', 0)} sampled runs; "
        f"photon ready: {photon}; dynamic SS: {dynamic}; approx. ohmic: {ohmic}.\n"
        "Ranges refer only to completed sampled currents; gray entries are unavailable, not physical failures."
    )


def _compact_range(item: Mapping[str, Any]) -> str:
    values = [float(value) for value in item.get("currents_uA", [])]
    if not values:
        return "none"
    if len(values) == 1:
        return rf"{values[0]:g} $\mu$A"
    return rf"{min(values):g}-{max(values):g} $\mu$A ({len(values)} samples)"


def _case_array(cases: Sequence[Mapping[str, Any]], field: str) -> np.ndarray:
    return np.asarray([float(case.get(field, np.nan)) for case in cases], dtype=float)


def _optional_bool(value: Any) -> bool | None:
    if value is None:
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return None


def _serializable_case(case: Mapping[str, Any]) -> dict[str, Any]:
    fields = (
        "run_name",
        "current_uA",
        "complete",
        "reason",
        "strict_stationarity_passes",
        "dynamic_stationarity_passes",
        "photon_ready",
        "photon_ready_source",
        "approximately_ohmic",
        "normal_voltage_ratio",
        "terminal_normal_voltage_ratio",
        "mean_delta_over_delta0",
        "normal_like_fraction_final",
        "accepted_steps",
        "rejected_steps",
        "rejected_over_accepted",
        "final_time_ps",
    )
    return {field: case.get(field) for field in fields if field in case}


__all__ = [
    "build_current_sweep_regime_summary",
    "plot_current_sweep_regime_summary",
    "write_current_sweep_regime_summary",
]
