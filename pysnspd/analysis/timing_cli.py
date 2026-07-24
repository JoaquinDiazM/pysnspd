"""Shared command-line options for photon timing post-processing."""

from __future__ import annotations

import argparse

from pysnspd.analysis.timing import DetectionCriteria, RecoveryCriteria


def add_timing_analysis_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--recovery-mode",
        choices=("electrical", "efficiency90", "state"),
        default="electrical",
    )
    parser.add_argument("--detection-threshold-uV", type=float, default=100.0)
    parser.add_argument(
        "--detection-polarity",
        choices=("positive", "negative", "auto"),
        default="positive",
    )
    parser.add_argument("--detection-baseline-window-ps", type=float, default=10.0)
    parser.add_argument("--detection-confirmation-ps", type=float, default=0.5)
    parser.add_argument("--detection-hysteresis-fraction", type=float, default=0.10)
    parser.add_argument("--peak-confirmation-ps", type=float, default=2.0)
    parser.add_argument("--post-peak-safety-ps", type=float, default=10.0)
    parser.add_argument("--recovery-hold-ps", type=float, default=10.0)
    parser.add_argument("--recovery-efficiency-fraction", type=float, default=0.90)
    parser.add_argument("--recovery-current-rel-tol", type=float, default=1.0e-2)
    parser.add_argument("--recovery-current-abs-uA", type=float, default=0.05)
    parser.add_argument("--recovery-voltage-rel-tol", type=float, default=1.0e-2)
    parser.add_argument("--recovery-voltage-abs-uV", type=float, default=10.0)
    parser.add_argument("--recovery-temperature-abs-K", type=float, default=0.05)
    parser.add_argument("--recovery-condensate-rel-tol", type=float, default=2.0e-2)
    parser.add_argument("--recovery-spatial-quantile", type=float, default=0.995)
    parser.add_argument("--recovery-spatial-max-guard-factor", type=float, default=4.0)


def timing_criteria_from_args(
    args: argparse.Namespace,
) -> tuple[DetectionCriteria, RecoveryCriteria]:
    return (
        DetectionCriteria(
            threshold_V=float(args.detection_threshold_uV) * 1.0e-6,
            polarity=str(args.detection_polarity),
            baseline_window_s=float(args.detection_baseline_window_ps) * 1.0e-12,
            confirmation_s=float(args.detection_confirmation_ps) * 1.0e-12,
            hysteresis_fraction=float(args.detection_hysteresis_fraction),
            peak_confirmation_s=float(args.peak_confirmation_ps) * 1.0e-12,
            post_peak_safety_s=float(args.post_peak_safety_ps) * 1.0e-12,
        ).validated(),
        RecoveryCriteria(
            mode=str(args.recovery_mode),
            hold_s=float(args.recovery_hold_ps) * 1.0e-12,
            efficiency_fraction=float(args.recovery_efficiency_fraction),
            current_relative_tolerance=float(args.recovery_current_rel_tol),
            current_absolute_tolerance_A=float(args.recovery_current_abs_uA) * 1.0e-6,
            voltage_relative_tolerance=float(args.recovery_voltage_rel_tol),
            voltage_absolute_tolerance_V=float(args.recovery_voltage_abs_uV) * 1.0e-6,
            temperature_absolute_tolerance_K=float(args.recovery_temperature_abs_K),
            condensate_relative_tolerance=float(args.recovery_condensate_rel_tol),
            spatial_quantile=float(args.recovery_spatial_quantile),
            spatial_max_guard_factor=float(args.recovery_spatial_max_guard_factor),
        ).validated(),
    )


__all__ = ["add_timing_analysis_arguments", "timing_criteria_from_args"]
