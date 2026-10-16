
'''Saved-Te reconstruction and U/G semantics; no simulation or full report.'''
import importlib.util
from pathlib import Path
import numpy as np
import pytest
from scipy.constants import Boltzmann as kB,hbar,elementary_charge as e
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("thermal_case_report",ROOT/"plot_pipelines/07_case_report.py")
report=importlib.util.module_from_spec(spec);spec.loader.exec_module(report)

def case(mode="electron"):
    Tc=8.65;Eu=kB*Tc;D=5e-5;Rs=608.;df=7e-9;ell=np.sqrt(hbar*D/(2*Eu))
    xy=np.array([[0.,0.],[2.,0.],[0.,1.]])*1e-9;edges=np.array([[0,1],[1,2],[2,0]])
    lens=np.linalg.norm(xy[edges[:,1]]-xy[edges[:,0]],axis=1);dual=np.array([.4,.7,.5])*1e-9
    mesh=dict(coordinates_m=xy,coordinates_bar=xy/ell,edges=edges,triangles=np.array([[0,1,2]]),
        edge_length_m=lens,dual_face_length_m=dual,conductance=dual/lens,
        area_weights=np.array([.3,.4,.3])*1e-18/ell**2,fixed_nodes=np.array([0,1,2]))
    m=dict(model_version="1.1.0-candidate",closure=dict(stiffness_temperature=mode),
        parameter_identity=dict(xi_cutoff=8.,xi_order=4),
        physicalparams=dict(energy_scale_J=Eu,N0_single_spin_per_J_m3=1/(2*e**2*D*Rs*df),
            film_thickness_m=df,bath_temperature_K=.9,diffusion_m2_s=D,sheet_resistance_ohm=Rs),
        geometry=dict(ell0_m=ell))
    return dict(m=m,mesh=mesh,dc=dict(left_nodes=np.array([0,2]),right_nodes=np.array([1])))

def test_candidate_gauge_reconstruction_uses_actual_saved_Te_and_free_G():
    c=case();f=report.reconstruct_functional(c)
    d=np.array([.8,1.2,.6])*np.exp(1j*np.array([.1,.7,-.4]));T=np.array([7.,12.,6.])
    p=f.fermi_dirac(d,T);ev=f.evaluate(d,p,temperature_K=T)
    exact=report.edge_supercurrent(f,d,T)
    np.testing.assert_allclose(exact,ev.edge_current_A,rtol=3e-15,atol=1e-30)
    cold=report.edge_supercurrent(f,d,np.full(3,.9))
    assert np.linalg.norm(exact-cold)>.2*np.linalg.norm(cold)
    with pytest.raises(ValueError,match="Saved Te_K"):
        report.edge_supercurrent(f,d)
    contract=report.model_contract(c["m"])
    assert "U_BCS+G-sum" in contract["stored_energy"]
    assert "fixed-Te" in contract["force_current"]
    assert "no capacity floor" in contract["admission"]
    assert "pysnspd/experimental/thermal_stiffness.py" in report.reconstruction_sources(c["m"])

def test_calibrated_D_bath_control_keeps_fixed_bath_reconstruction():
    c=case("bath");c["m"]["physicalparams"]["diffusion_m2_s"]=1.4e-4
    f=report.reconstruct_functional(c)
    d=np.array([.8,1.2,.6])*np.exp(1j*np.array([.1,.7,-.4]))
    np.testing.assert_array_equal(report.edge_supercurrent(f,d),report.edge_supercurrent(f,d,np.full(3,20.)))
    assert report.model_contract(c["m"])["bath_control"]
    assert not report.model_contract(c["m"])["candidate"]
    assert report.reconstruction_sources(c["m"])==report.RECON_SOURCES

def test_saved_snapshot_Te_recovers_Vdev_and_current_maps():
    c=case();f=report.reconstruct_functional(c);d=np.array([.8,1.2,.6])*np.exp(1j*np.array([.1,.7,-.4]))
    T=np.array([7.,12.,6.]);Is=2e-6
    ports=report.EffectiveNeutralPorts(f.graph,608.,c["dc"]["left_nodes"],c["dc"]["right_nodes"])
    el=ports.solve(f.evaluate(d,f.fermi_dirac(d,T),temperature_K=T).edge_current_A,Is)
    c["s"]=dict(time_s=np.array([0.]),delta_J=np.array([d*f.Eu_J]),Te_K=np.array([T]))
    c["z"]=dict(time_s=np.array([0.]),Is_A=np.array([Is]),Vdev_V=np.array([el.Vdev_V]))
    fields,checks=report.currents(c,np.array([0]))
    assert checks[0]["Vdev_error_V"]<1e-16
    assert checks[0]["model_contract"]["stiffness_temperature"]=="saved Te_K"
    assert np.all(np.isfinite(fields[0]["total"]))
    c["s"]["Te_K"][:]=.9
    with pytest.raises(ArithmeticError,match="Vdev differs"):
        report.currents(c,np.array([0]))

def test_thermal_identity_cannot_be_inferred_from_an_ambiguous_version():
    c=case();c["m"].pop("model_version")
    with pytest.raises(ValueError,match="Explicit model_version"):
        report.reconstruct_functional(c)

@pytest.mark.parametrize("version",["1.1.0-candidate","1.1.1-candidate","1.1.2-candidate"])
def test_candidate_summary_equations_render_and_state_internal_energy_contract(version):
    class Capture:
        def __init__(self):
            self.c=case();self.c["m"]["model_version"]=version;self.caption=None;self.fig=None
        def page(self,title):
            self.title=title;self.fig=report.plt.figure(figsize=(11.69,8.27));return self.fig
        def save(self,fig,title,caption):
            self.caption=caption;fig.canvas.draw()
    r=Capture()
    try:
        report.model_summary(r)
        assert "candidate" in r.title
        assert "stored internal energy differs from free G" in r.caption
        assert "Ceff>0" in r.caption
        # Equations stay within the landscape page rather than being clipped.
        renderer=r.fig.canvas.get_renderer()
        for text in r.fig.texts:
            box=text.get_window_extent(renderer=renderer)
            assert box.x1<=r.fig.bbox.x1+1
    finally:report.plt.close(r.fig)

@pytest.mark.parametrize("version",["1.1.1-candidate","1.1.2-candidate"])
def test_current_reports_use_Pauli_heat_contract(version):
    m={"model_version":version,"closure":{"stiffness_temperature":"electron"}}
    assert report.thermal_candidate(m)
    assert report.model_contract(m)["pauli_signed_heat"]
    old={"model_version":"1.1.0-candidate","closure":{"stiffness_temperature":"electron"}}
    assert not report.model_contract(old)["pauli_signed_heat"]
