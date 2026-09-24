"""Build the auditable, user-executed pre-photon DC campaign; no simulation."""
from pathlib import Path
import json
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from sandbox.stage4_core.coupled_campaign import atomic_json,sha
from sandbox.stage5_prephoton.provenance import source_sha

OUT=ROOT/'docs/implementation/stage5/prephoton_dc_20260924'
BASE='docs/implementation/stage4/practical_time_review_20260924/dual_mesh/resampled/mesh.npz'

def reference(mesh,current):
    return dict(schema='pysnspd.stage5.intrinsic_dc_reference.v1',mesh=mesh,mesh_sha256=sha(ROOT/mesh),
        T_K=.9,Tc_K=8.65,diffusion_m2_s=5e-5,sheet_resistance_ohm=608.,width_m=80e-9,
        current_A=current,matsubara_count=256,active_length_m=5e-6,added_inductance_H=96e-9,
        spectral_tolerance=1e-8,max_newton_iterations=100,max_outer_iterations=1000,
        gap_rms_relative_tolerance=2e-6,bulk_gap_relative_limit=.01,current_relative_limit=.02)

def main():
    directory=OUT/'plans';directory.mkdir(parents=True,exist_ok=True)
    definitions=[('dc155_l160_h5',BASE,15.5e-6),('dc215_l160_h5',BASE,21.5e-6),
        ('dc215_l160_h35',str((OUT/'meshes/l160_h35/mesh.npz').relative_to(ROOT)).replace('\\','/'),21.5e-6),
        ('dc215_l320_h5',str((OUT/'meshes/l320_h5/mesh.npz').relative_to(ROOT)).replace('\\','/'),21.5e-6)]
    references=[]
    for name,mesh,current in definitions:
        path=directory/(name+'.json');atomic_json(path,reference(mesh,current))
        references.append(dict(id=name,plan=path.relative_to(ROOT).as_posix(),sha256=sha(path),
            estimated_seconds=1800 if name.endswith(('h35','l320_h5')) else 900))
    template=dict(schema='pysnspd.stage5.prephoton_dc_hold.v1',T_K=.9,Tc_K=8.65,
        tau_ee_Tc_ps=6.,tau_ep_Tc_ps=24.7,delta0_over_kBTc=1.764,
        sheet_resistance_ohm=608.,duration_ps=10.,spectral_tolerance=1e-8,
        maximum_newton_iterations=100,gap_drift_relative_limit=.01,gap_ripple_relative_limit=.01,
        current_cut_relative_limit=.02,current_absolute_floor_A=1e-12,
        maximum_relative_displacement=.1,observation_cadence=100,
        dark_readout_fraction_limit=.001,energy_balance_relative_limit=.01,
        circuit=dict(R_bias_ohm=1e4,L_bias_H=1e-6,R_load_ohm=50.,C_couple_F=100e-12))
    holds=[dict(id=item['id'],reference=item['id'],step_multiplier=1.,
        estimated_seconds=28800 if item['id'].endswith(('h35','l320_h5')) else 14400) for item in references]
    holds.append(dict(id='dc215_l160_h5_half_step',reference='dc215_l160_h5',step_multiplier=.5,estimated_seconds=28800))
    files=['pysnspd/experimental/bulk_current_reference.py','pysnspd/experimental/dc_current_ports.py',
        'pysnspd/experimental/thermal_spatial_usadel.py','pysnspd/experimental/thermal_stable_newton.py',
        'pysnspd/experimental/thermal_weak_response.py','pysnspd/experimental/heredado_kwt_bridge.py',
        'pysnspd/experimental/stable_kwt_euler.py',
        'pysnspd/experimental/thermal_snapshot.py','pysnspd/experimental/energy_catalog.py',
        'pysnspd/experimental/electrical_ports.py','pysnspd/solver/core.py',
        'sandbox/stage4_core/coupled_campaign.py','sandbox/stage4_core/parallel_runtime.py',
        'sandbox/stage4_core/dual_kwt_time.py','sandbox/stage4_core/biased_strip_reference.py',
        'sandbox/stage4_core/run_self_consistent_core.py','sandbox/stage4_core/biased_coupled_response.py',
        'sandbox/stage4_core/coupled_response.py',
        'sandbox/stage5_prephoton/prepare_dc_reference.py','sandbox/stage5_prephoton/dc_hold.py',
        'sandbox/stage5_prephoton/provenance.py',
        'sandbox/stage5_prephoton/run_prephoton_dc.py','sandbox/stage5_prephoton/analyze_prephoton_dc.py']
    plan=dict(schema='pysnspd.stage5.prephoton_dc_campaign.v1',references=references,holds=holds,
        hold_template=template,maximum_parallel_cases=2,maximum_workers_per_case=13,
        euler_safety=.4,maximum_dt_ps=.0001,
        scope='Current-carrying bulk cuts, DC only, fixed external inductance, no photon; bath equilibrium preserved, not nonlinear heat validation.',
        source_hash_normalization='CRLF-to-LF',source_sha256={p:source_sha(ROOT/p) for p in files},
        comparisons=dict(mesh=['dc215_l160_h5','dc215_l160_h35'],length=['dc215_l160_h5','dc215_l320_h5'],
            time=['dc215_l160_h5','dc215_l160_h5_half_step']),
        acceptance=dict(gap_comparison_relative=.01,current_comparison_relative=.02,normalized_dark_readout_difference=.001),
        estimate_note='Planning estimates only. Per-case ETA is replaced by measured progress; no guaranteed parallel speedup.')
    atomic_json(OUT/'campaign.json',plan)
    print(json.dumps(dict(plan=str(OUT/'campaign.json'),references=len(references),holds=len(holds))))

if __name__=='__main__':main()
