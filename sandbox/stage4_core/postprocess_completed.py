"""Algebraic boundary-response review of immutable stage4A fields; no spectra."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from controls import catalogue, material_reference, scales
from pysnspd.experimental.cell_closures import KWTMobility
from pysnspd.experimental.energy_catalog import K_B_J_K
from pysnspd.experimental.spatial_dynamics import kwt_spatial_response


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-root', type=Path, default=ROOT/'docs/implementation/stage4/review_20260923/raw')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    plan = json.loads((ROOT/'docs/implementation/stage4/start_20260923/campaign_plan.json').read_text())
    cat = catalogue(plan, 'reference')  # Construction only: zero energy-kernel queries.
    m = material_reference(plan)
    data = {}
    for path in sorted((args.raw_root/'stage4A_spatial_20260923').glob('*/fields.npz')):
        with np.load(path) as a:
            z, g, xy = a['delta'], a['cartesian_force'], a['coordinates_m']
            mass = a['quadrature_volume_m3']/(m['width_m']*m['thickness_m']*m['ell0_m'])
            temperature = a['temperature_equivalent_K']*K_B_J_K/cat.vacuum.delta0_J
            # Coordinates are stored exactly; all outer rows, including corners once.
            boundary = np.any((xy == xy.min(axis=0)) | (xy == xy.max(axis=0)), axis=1)
            load = np.zeros_like(g)
            load[boundary] = g[boundary]
            record = dict(source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                          nodes=len(z), constrained_nodes=int(boundary.sum()), mobility={})
            for name, pair in plan['mobility_pairs_ps'].items():
                mobility = KWTMobility(scales(cat, plan), *pair)
                free = kwt_spatial_response(z, g, mass, mobility, temperature, np.zeros(len(z)))
                fixed = kwt_spatial_response(z, g, mass, mobility, temperature, np.zeros(len(z)),
                                             boundary_load_bar=load)
                record['mobility'][name] = dict(
                    original_heat_density_max_reproduction_error=float(np.max(abs(free.heat_density_bar-a['heat_'+name]))),
                    original_total_heat=free.condensate_heat_rate_bar,
                    original_boundary_heat=float(np.dot(mass[boundary], free.heat_density_bar[boundary])),
                    original_interior_heat=float(np.dot(mass[~boundary], free.heat_density_bar[~boundary])),
                    constrained_total_heat=fixed.condensate_heat_rate_bar,
                    constrained_field_work=fixed.field_energy_rate_bar,
                    constraint_work=fixed.boundary_work_rate_bar,
                    identity_residual=fixed.identity_residual_bar,
                    maximum_boundary_velocity=float(np.max(abs(fixed.field_velocity_bar[boundary]))),
                    interior_velocity_max_change=float(np.max(abs(fixed.field_velocity_bar[~boundary]-free.field_velocity_bar[~boundary]))),
                    reaction_integrated_norm=float(np.linalg.norm(load)),
                    reaction_rule='Applied conjugate load equals G on prescribed boundary nodes; zero elsewhere')
            data[path.parent.name] = record
    result = dict(schema='pysnspd.stage4A.boundary_review.v1', cases=data,
        scope='New algebraic response on archived fields. Four sides clamped, phi=0; not a reservoir/device boundary model or time trajectory.',
        reason='For diagonal nodal KWT at phi=0, load=G on fixed nodes cancels their material velocity without deleting energy-work terms.',
        energy_kernel_queries=cat.calls, stage4_complete=False, photon=False)
    assert len(data) == 6 and cat.calls == 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf8')
    print(json.dumps(dict(cases=len(data), energy_kernel_queries=cat.calls, output=str(args.output))))


if __name__ == '__main__':
    main()
