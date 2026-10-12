import copy,importlib.util,sys
from pathlib import Path
import pytest
BASE=Path(__file__).parent
if not (BASE/"width_guard_window_candidate.py").exists():BASE=BASE.parent/"sandbox/stage5_effective_validation"
if str(BASE) not in sys.path:sys.path.insert(0,str(BASE))
s=importlib.util.spec_from_file_location("width_guard_window_candidate",BASE/"width_guard_window_candidate.py")
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def histories():
    ref=[dict(start_time_ps=0.,end_time_ps=12.4,local_energy_tolerance=None),
        dict(start_time_ps=12.4,end_time_ps=22.435,local_energy_tolerance=1e-5),
        dict(start_time_ps=22.435,end_time_ps=None,local_energy_tolerance=2.5e-6)]
    fresh=[dict(start_time_ps=0.,end_time_ps=None,local_energy_tolerance=1e-5)]
    return ref,fresh
def test_guard_is_compared_in_original_window_not_from_final_parameter():
    r,c=histories();law=m.window_guard_law(r,c,[0.,20.])
    assert law["passed"] and law["numerical_admission_axis"] is False
    assert law["policy_intervals"][0]["earlier_legacy_nominal_guard_disabled"] is True
    assert law["policy_intervals"][1]["nominal_guard"]==1e-5
@pytest.mark.parametrize("interval",[[0.,200.],[10.,20.],[0.,22.435]])
def test_only_explicit_original20ps_window_supported(interval):
    with pytest.raises(ValueError):m.window_guard_law(*histories(),interval)
def test_reference_transition_inside_width_window_must_be_reviewed():
    r,c=histories();r[1]["end_time_ps"]=15.;r[2]["start_time_ps"]=15.
    with pytest.raises(ValueError):m.window_guard_law(r,c,[0.,20.])
def test_width_policy_cannot_be_looser_or_replaced_by_later_guard():
    r,c=histories()
    for value in (2e-5,2.5e-6,None):
        c[0]["local_energy_tolerance"]=value
        with pytest.raises(ValueError):m.window_guard_law(r,c,[0.,20.])
def case(label="narrow"):
    sigma=7.549834435;factor=m.RATIOS[label]
    p=dict(sigma_nm=sigma,local_energy_tolerance=2.5e-6,dt_ps=.02,step_tolerance=2e-7,tau_ep_Tc_ps=24.7)
    source=dict(wavelength_nm=775.,photon_energy_J=2.,retained_fraction=.667,retained_energy_J=1.,
        center_m=[1.,0.],time_origin="postcascade",initialization="same",sigma_nm=sigma)
    ref=dict(manifest=dict(parameter_identity=p,source_pins_sha256={str(i):str(i) for i in range(13)},
        physicalparams=dict(Tc_K=8.65),source=source))
    candidate=copy.deepcopy(ref);candidate["manifest"]["parameter_identity"].update(sigma_nm=sigma*factor,local_energy_tolerance=1e-5)
    candidate["manifest"]["source"]["sigma_nm"]=sigma*factor
    return ref,candidate
@pytest.mark.parametrize("label",["narrow","wide"])
def test_only_declared_physical_sigma_factor_with_guard_window_exception(label):
    m.validate_physical_parameters(*case(label),label,{"threads","output_dir"})
@pytest.mark.parametrize("fault",["tau","dt","solver","source","pins","circuit","sigma","guard"])
def test_no_general_parameter_ignore_for_width_cases(fault):
    r,c=case();cm=c["manifest"]
    if fault=="tau":cm["parameter_identity"]["tau_ep_Tc_ps"]=5.
    if fault=="dt":cm["parameter_identity"]["dt_ps"]=.1
    if fault=="solver":cm["parameter_identity"]["step_tolerance"]=2e-6
    if fault=="source":cm["source"]["retained_energy_J"]=2.
    if fault=="pins":cm["source_pins_sha256"]["0"]="changed"
    if fault=="circuit":cm["physicalparams"]["Tc_K"]=9.
    if fault=="sigma":cm["parameter_identity"]["sigma_nm"]*=2.
    if fault=="guard":cm["parameter_identity"]["local_energy_tolerance"]=2e-5
    with pytest.raises(ValueError):m.validate_physical_parameters(r,c,"narrow",{"threads","output_dir"})
