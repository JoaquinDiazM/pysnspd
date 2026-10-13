"""Physical reconstruction and per-case report semantics; never evolve a model."""
import importlib.util
from pathlib import Path
import json
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("case_report",ROOT/"plot_pipelines/07_case_report.py")
report=importlib.util.module_from_spec(spec);spec.loader.exec_module(report)


def small_mesh():
    xy=np.array([[0.,0.],[2.,0.],[0.,1.]])*1e-9
    edges=np.array([[0,1],[1,2],[2,0]])
    return dict(coordinates_m=xy,edges=edges,triangles=np.array([[0,1,2]]),
        edge_length_m=np.linalg.norm(xy[edges[:,1]]-xy[edges[:,0]],axis=1),
        dual_face_length_m=np.array([.4,.7,.5])*1e-9)


def test_exact_gauge_current_is_independent_of_unsaved_occupations():
    from pysnspd.experimental.thermal_spatial_usadel import ThermalGraph
    from pysnspd.experimental.effective_bcs_functional import EffectiveBCSFunctional,bcs_count_quadrature
    mesh=small_mesh()
    graph=ThermalGraph(np.ones(3),mesh["edges"],np.ones(3),mesh["coordinates_m"]/1e-9,np.array([],int))
    xi,w=bcs_count_quadrature(8.,4)
    f=EffectiveBCSFunctional(graph,energy_scale_J=1.194261385e-22,N0_per_J_m3=9.1533e46,
        ell0_m=4.69849e-9,thickness_m=7e-9,count_nodes=xi,count_weights=w,
        gradient_law="bcs_radial",bath_temperature_K=.9,diffusion_m2_s=5e-5)
    d=np.array([.8,1.2,.6])*np.exp(1j*np.array([.1,.7,-.4]))
    p=np.random.default_rng(4).uniform(0,.4,(3,len(xi)))
    exact=report.edge_supercurrent(f,d)
    np.testing.assert_allclose(exact,f.evaluate(d,p).edge_current_A,rtol=2e-15,atol=1e-30)
    np.testing.assert_allclose(exact,f.evaluate(d,np.zeros_like(p)).edge_current_A,rtol=2e-15,atol=1e-30)
    np.testing.assert_allclose(exact,report.edge_supercurrent(f,d*np.exp(1.234j)),rtol=2e-14,atol=1e-30)


def test_density_reconstruction_recovers_uniform_vector_and_orientation():
    mesh=small_mesh();target=np.array([3e10,-2e10]);df=7e-9
    unit=(mesh["coordinates_m"][mesh["edges"][:,1]]-mesh["coordinates_m"][mesh["edges"][:,0]])/mesh["edge_length_m"][:,None]
    flux=(unit@target)*mesh["dual_face_length_m"]*df
    vectors,check=report.triangle_vectors(mesh,flux,df)
    np.testing.assert_allclose(vectors,[target],rtol=1e-14)
    assert check["relative_projection_RMS"]<1e-14
    reversed_mesh={**mesh,"edges":mesh["edges"][:,::-1]}
    reversed_vectors,_=report.triangle_vectors(reversed_mesh,-flux,df)
    np.testing.assert_allclose(vectors,reversed_vectors,rtol=1e-14)
    other=np.array([-1e10,4e10]);flux2=(unit@other)*mesh["dual_face_length_m"]*df
    summed,_=report.triangle_vectors(mesh,flux+flux2,df)
    np.testing.assert_allclose(summed,[target+other],rtol=1e-14)


def test_zero_dual_face_is_rejected():
    mesh=small_mesh();mesh["dual_face_length_m"][0]=0
    with pytest.raises(ValueError,match="Positive dual"):report.triangle_vectors(mesh,np.ones(3),7e-9)


def test_hotbelt_requires_actual_connected_transverse_path():
    xy=np.array([[0,-1],[0,0],[0,1],[1,-1],[1,0],[1,1]],float)
    edges=np.array([[0,1],[1,2],[3,4],[4,5],[0,3],[1,4],[2,5]])
    fields=np.array([[.1,1,.1,1,1,1],[.1,.1,.1,1,1,1]],complex)
    assert report.hotbelts(fields,np.ones(6),edges,xy).tolist()==[False,True]


def metric_case(provisional=False,detected=False,censored=False):
    t=np.linspace(0,1e-12,11);voltage=t/1e-12*.2e-3 if detected else t/1e-12*1e-6
    z={"time_s":t,"Vout_V":voltage,"Ib_A":np.full(11,20e-6),"Is_A":20e-6-voltage/50}
    obs=dict(trigger_confirmed=detected,censored=censored,baseline_Vout_V=0.,
        first_crossing_time_s=.5e-12 if detected else None,
        confirmation_time_s=.8e-12 if detected else None,config={"event_time_s":0.})
    return dict(z=z,m={"physicalparams":{"circuit":{"R_load_ohm":50.}}},obs=obs,
        s={"time_s":t},provisional=provisional)


def test_confirmed_latency_and_physical_local_slope():
    result=report.metrics(metric_case(detected=True))
    assert result["latency_ps"]==pytest.approx(.5)
    assert result["confirmation_ps"]==pytest.approx(.8)
    assert result["local_slope_mV_ps"]==pytest.approx(.2)
    assert result["readout_identity_error_V"]<1e-18


@pytest.mark.parametrize("partial,censor,expected",[(False,True,True),(True,True,False),(True,False,False)])
def test_provisional_is_never_terminal_censoring(partial,censor,expected):
    result=report.metrics(metric_case(provisional=partial,censored=censor))
    assert result["censored"] is expected
    assert result["latency_ps"] is None
    assert result["local_slope_mV_ps"] is None


def test_broken_readout_identity_fails():
    case=metric_case();case["z"]["Vout_V"]*=2
    with pytest.raises(ValueError,match="identity"):report.metrics(case)



def test_snapshot_selection_always_has_six_ordered_slots():
    t=np.array([0,.1,1,10,15,20,30,31,32])*1e-12
    obs=dict(trigger_confirmed=True,first_crossing_time_s=30.4e-12,confirmation_time_s=30.9e-12)
    hb=np.array([False]*4+[True]*5)
    auto=report.select_snapshots(t,obs,hb)
    assert auto.tolist()==[0,2,4,6,7,8]
    explicit=report.select_snapshots(t,{},hb,[14.9,.16])
    assert explicit.tolist()==[4,1,1,1,1,1]
    assert report.snapshot_targets(t,{},hb,[14.9,.16]).tolist()==[14.9,.16,.16,.16,.16,.16]


def test_snapshot_extras_are_ignored_without_reordering_or_endpoints():
    times=np.arange(8)*1e-12
    requested=[5,4,3,2,1,.2,7,float("nan")]
    indices=report.select_snapshots(times,{},np.zeros(8,bool),requested)
    assert indices.tolist()==[5,4,3,2,1,0]
    assert report.select_snapshots(times,{},np.zeros(8,bool),[4]*6).tolist()==[4]*6


def test_snapshot_nearest_endpoint_and_sparse_acquisition_receipt():
    times=np.array([0,1,2])*1e-12
    targets=report.snapshot_targets(times,{},np.zeros(3,bool),[-2,.6,100])
    ids=report.select_snapshots(times,{},np.zeros(3,bool),[-2,.6,100])
    assert ids.tolist()==[0,1,2,2,2,2]
    receipt=report.snapshot_selection(times,targets,ids)
    assert receipt["slots"][0]["outside_saved_window"] is True
    assert receipt["slots"][1]["saved_ps"]==1
    assert receipt["slots"][3]["repeated_saved_state"] is True
    assert receipt["future_visualization_queue"]
    assert report.select_snapshots(np.array([0.]),{},np.zeros(1,bool),[5]).tolist()==[0]*6


@pytest.mark.parametrize("requested",[[],[float("nan")],[float("inf")]])
def test_nonfinite_or_empty_used_snapshot_requests_fail(requested):
    with pytest.raises(ValueError,match="finite requested"):
        report.select_snapshots(np.array([0,1e-12]),{},np.zeros(2,bool),requested)


def vortex_triangle():
    xy=np.array([[-1,-1],[1,-1],[0,1]],float)*1e-9
    delta=(xy[:,0]+1j*xy[:,1])/1e-9*1e-22
    return xy,delta,np.array([[0,1,2]])


def test_phase_vortex_sign_core_and_global_phase_invariance():
    xy,delta,tri=vortex_triangle();dc=np.full(3,1e-22)
    vortex=report.phase_vorticity(delta,dc,xy,tri)
    assert vortex["positive"]==1 and vortex["negative"]==0
    assert vortex["valid"].tolist()==[True]
    np.testing.assert_allclose(vortex["core_points_m"],[[0,0]],atol=1e-24)
    anti=report.phase_vorticity(np.conj(delta),dc,xy,tri)
    assert anti["negative"]==1 and anti["positive"]==0
    rotated=report.phase_vorticity(delta*np.exp(2.3j),dc,xy,tri[:,::-1])
    assert rotated["positive"]==1
    np.testing.assert_allclose(rotated["core_points_m"],vortex["core_points_m"],atol=1e-24)


def test_phase_gradient_is_not_vorticity_and_zero_contours_are_undefined():
    xy,delta,tri=vortex_triangle();dc=np.full(3,1e-22)
    smooth=report.phase_vorticity(np.exp(1j*xy[:,0]/1e-8)*1e-22,dc,xy,tri)
    assert smooth["winding"].tolist()==[0] and smooth["valid"].all()
    nodal_zero=delta.copy();nodal_zero[0]=0
    ambiguous=report.phase_vorticity(nodal_zero,dc,xy,tri)
    assert ambiguous["indeterminate_cells"]==1
    assert not len(ambiguous["core_points_m"])
    suppressed=report.phase_vorticity(delta*1e-10,dc,xy,tri)
    assert suppressed["indeterminate_cells"]==1


def test_vortex_antivortex_pair_preserves_signed_counts():
    from scipy.spatial import Delaunay
    x,y=np.meshgrid(np.linspace(-2,2,17),np.linspace(-2,2,17))
    xy=np.column_stack((x.ravel(),y.ravel()))*1e-9
    z=xy[:,0]/1e-9+1j*xy[:,1]/1e-9
    delta=(z-(-.63+.17j))*np.conj(z-(.61-.23j))*1e-22
    result=report.phase_vorticity(delta,np.full(len(xy),1e-22),xy,Delaunay(xy).simplices)
    assert (result["positive"],result["negative"])==(1,1)
    assert result["core_charge"].sum()==0


def test_quantity_families_have_distinct_palettes():
    assert report.PALETTES["energy"]!=report.PALETTES["signed_current"]
    assert report.PALETTES["energy"]!=report.PALETTES["total_current"]
    assert len(set(report.PALETTES.values()))==len(report.PALETTES)


def test_decompression_budget_uses_uncompressed_sizes(tmp_path):
    p=tmp_path/"zeros.npz";np.savez_compressed(p,field=np.zeros(1000000))
    assert p.stat().st_size<100000
    with pytest.raises(MemoryError):
        report.memory_check([(p,["field"])],513*1024**2)


def test_same_basename_in_different_runs_has_unique_output_identity(tmp_path):
    assert report.case_id(tmp_path/"a"/"run")!=report.case_id(tmp_path/"b"/"run")


def fixture_run(tmp_path,complete):
    from pysnspd.experimental.trigger_window import TriggerWindow,TriggerWindowConfig
    run=tmp_path/"run";(run/"dc").mkdir(parents=True)
    xy=np.array([[0,-1],[1,-1],[0,1],[1,1]],float)*1e-9
    edges=np.array([[0,1],[1,3],[3,2],[2,0],[0,3]])
    triangles=np.array([[0,1,3],[0,3,2]])
    meshfile=tmp_path/"mesh.npz"
    np.savez(meshfile,coordinates_m=xy,coordinates_bar=xy/1e-9,triangles=triangles,edges=edges,
        area_weights=np.ones(4),conductance=np.ones(5),dual_area_m2=np.ones(4)*1e-18,
        edge_length_m=np.ones(5)*1e-9,dual_face_length_m=np.ones(5)*1e-9,fixed_nodes=np.arange(4))
    dcfile=run/"dc/effective_dc.npz"
    np.savez(dcfile,delta_bar=np.ones(4),coordinates_m=xy,left_nodes=np.array([0,2]),right_nodes=np.array([1,3]))
    (run/"dc/effective_dc_receipt.json").write_text(json.dumps({"fields_sha256":report.sha256(dcfile)}))
    t=np.array([0,1e-12]);trace={key:np.zeros(2) for key in report.TRACE};trace["time_s"]=t
    trace["Ib_A"][:]=20e-6;trace["Is_A"][:]=20e-6
    maps={key:np.ones((2,4)) for key in report.MAPS[1:]};maps["delta_J"]=np.ones((2,4),complex)*1e-22
    maps["time_s"]=t
    pins={p:report.sha256(ROOT/p) for p in report.RECON_SOURCES}
    identity={"fixture":True}
    controller=TriggerWindow(TriggerWindowConfig(0.,0.,"fixture DC",1e-4,"positive",.5e-12,.05,5e-12,200e-12))
    controller.observe_accepted(0.,0.);controller.observe_accepted(1e-12,0.)
    manifest=dict(schema="pysnspd.effective_photon_run.v1",closure={"microscopic_exact_Usadel":False},
        source_pins_sha256=pins,parameter_identity=identity,source={},run_complete=complete,
        status="TIME_CENSORED_NO_TRIGGER" if complete else "RUNNING",
        observation={"baseline_Vout_V":0.},geometry={"source_path":str(meshfile),"mesh_sha256":report.sha256(meshfile)})
    np.savez(run/"trajectory.npz",**trace)
    np.savez(run/"snapshots.npz",**maps,coordinates_m=xy,triangles=triangles,edges=edges)
    manifest["output_sha256"]={p:report.sha256(run/p) for p in ("trajectory.npz","snapshots.npz","dc/effective_dc.npz")}
    (run/"manifest.json").write_text(json.dumps(manifest))
    meta=dict(parameter_identity=identity,source_pins_sha256=pins,time_ps=1.,source={},controller=controller.checkpoint())
    np.savez(run/"checkpoint_latest.npz",metadata_json=np.array(json.dumps(meta)),
        **{"trajectory__"+k:v for k,v in trace.items()},**{"snapshots__"+k:v for k,v in maps.items()})
    return run


def test_incomplete_reads_one_checkpoint_and_preserves_running_input(tmp_path):
    run=fixture_run(tmp_path,False);build=tmp_path/"build";build.mkdir()
    before=report.sha256(run/"checkpoint_latest.npz")
    case=report.load_case(run,build,True,2**32)
    assert case["provisional"] is True
    assert case["obs"]["detected"] is False
    assert case["obs"]["censored"] is False
    assert case["z"]["time_s"][-1]==1e-12
    assert report.sha256(build/"checkpoint_capture.npz")==before
    assert report.sha256(run/"checkpoint_latest.npz")==before


def test_incomplete_requires_explicit_opt_in(tmp_path):
    run=fixture_run(tmp_path,False);build=tmp_path/"build";build.mkdir()
    with pytest.raises(ValueError,match="Incomplete"):report.load_case(run,build,False,2**32)
    assert not(build/"checkpoint_capture.npz").exists()


def test_completed_input_corruption_is_rejected(tmp_path):
    run=fixture_run(tmp_path,True);build=tmp_path/"build";build.mkdir()
    (run/"trajectory.npz").write_bytes(b"changed")
    with pytest.raises(ValueError,match="checksum"):report.load_case(run,build,False,2**32)


def test_physical_source_change_is_not_silently_reconstructed(tmp_path):
    run=fixture_run(tmp_path,True);build=tmp_path/"build";build.mkdir()
    p=run/"manifest.json";m=json.loads(p.read_text());m["source_pins_sha256"][report.RECON_SOURCES[0]]="wrong"
    p.write_text(json.dumps(m))
    with pytest.raises(ValueError,match="source differs"):report.load_case(run,build,False,2**32)


def test_fixed_a4_page_ignores_thesis_tight_bbox(monkeypatch):
    import matplotlib.pyplot as plt
    captured=[]
    class PDF:
        def savefig(self,figure,**kwargs):
            captured.append((figure.get_size_inches().tolist(),plt.rcParams["savefig.bbox"]))
    r=report.Report.__new__(report.Report)
    r.c={"metrics":{"final_ps":1.}};r.pdf=PDF();r.pages=[];r.dpi=160
    with plt.rc_context({"savefig.bbox":"tight"}):
        r.save(plt.figure(figsize=(11.69,8.27)),"test","caption")
        assert plt.rcParams["savefig.bbox"]=="tight"
    assert captured==[([11.69,8.27],None)]


def test_physical_voltage_unit_selection():
    assert report.voltage_units(np.array([.02]))==(1e3,"mV")
    assert report.voltage_units(np.array([36e-6]))==(1e6,"µV")
    assert report.voltage_units(np.array([9e-9]))==(1e9,"nV")


def test_worker_initialization_allocates_distinct_cores(monkeypatch):
    from multiprocessing import Value
    counter=Value("i",0);calls=[]
    monkeypatch.setattr(report.os,"sched_setaffinity",lambda pid,cpus:calls.append(cpus))
    report.initialize_worker([12,11],counter)
    report.initialize_worker([12,11],counter)
    assert calls==[{12},{11}]


def test_photon_free_hold_is_not_labeled_as_photon(tmp_path):
    run=tmp_path/"dc";run.mkdir()
    (run/"manifest.json").write_text(json.dumps(dict(schema="pysnspd.effective_photon_run.v1",closure=dict(microscopic_exact_Usadel=False),parameter_identity=dict(dc_only=True))))
    with pytest.raises(ValueError,match="Photon case required"):
        report.load_case(run,tmp_path,False,2**30)


def test_snapshot_layout_keeps_six_axes_with_larger_labels():
    import matplotlib.pyplot as plt
    fig=plt.figure(figsize=(11.69,8.27))
    axes=report.snapshot_axes(fig)
    assert len(axes)==6
    assert len({round(ax.get_position().y0,2) for ax in axes})==2
    plt.close(fig)


def test_contour_crossing_zero_is_not_a_resolved_vortex():
    xy,delta,tri=vortex_triangle()
    delta=np.array([-1,1,1+1j])*1e-22
    d=report.phase_vorticity(delta,np.full(3,1e-22),xy,tri)
    assert d["indeterminate_cells"]==1
    assert d["positive"]==d["negative"]==0
