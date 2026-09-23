"""Assemble a review ledger; ranges are scenarios/proposals, not fitted confidence intervals."""
from pathlib import Path
import json,hashlib,math
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'docs/implementation/stage3_5/research_20260923'
BASE=ROOT/'docs/implementation/stage3_5/parameter_inventory.json'

def main():
    inventory=json.loads(BASE.read_text(encoding='utf8'))
    domains={
      'material':('material_kinetics.md','MATERIAL_EVIDENCE_CONDITIONAL'),
      'bath':('material_kinetics.md','REFERENCE_SCENARIO'),
      'geometry':('geometry_numerics.md','CONDITIONAL_DESIGN'),
      'photon':('photon_preparation.md','TRANSFER_NOT_MEASURED'),
      'electronic_state':('geometry_numerics.md','DYNAMIC_STATE_NOT_INDEPENDENT_INPUT'),
      'condensate':('material_kinetics.md','EFFECTIVE_MODEL_CONSTRAINT'),
      'kinetic':('material_kinetics.md','MATERIAL_CLOSURE_NOT_CALIBRATED'),
      'boundaries':('geometry_numerics.md','COUPLED_BOUNDARY_REQUIREMENT'),
      'circuit':('circuit_observables.md','MEMORY_CIRCUIT_SCENARIO'),
      'catalogue':('geometry_numerics.md','NUMERICAL_COVERAGE_NOT_PHYSICAL_RANGE'),
      'spatial_numerics':('geometry_numerics.md','NUMERICAL_PROPOSAL_REQUIRES_TRAJECTORY'),
      'time_numerics':('geometry_numerics.md','NUMERICAL_PROPOSAL_REQUIRES_TRAJECTORY'),
      'observation':('circuit_observables.md','OBSERVABLE_DEFINITION_OR_MISSING_MEASUREMENT'),
      'provenance':('README.md','TRACEABILITY_REQUIREMENT'),
      'experiment_reference':('circuit_observables.md','REFERENCE_MUST_RETAIN_SAMPLE_IDENTITY')}
    specific={}
    def setrow(ids,status,value,note):
        for i in ids.split():specific[i]=(status,value,note)
    e=1.602176634e-19;kB=1.380649e-23;hbar=1.054571817e-34
    D=5e-5;rho=608*7e-9;sigma=1/rho;N0=sigma/(2*e*e*D)
    setrow('Tc','REFERENCE_SCENARIO',{'value':8.65,'units':'K'},'K20 reference film; no measurement uncertainty supplied.')
    setrow('D','USER_SELECTED_ADJUSTED_REFERENCE',{'value':D,'units':'m2/s','legacy_comparison':1.581e-4},'Transport spread parameter. K20 model input, not independently measured on the selected nanowire.')
    setrow('sheet_resistance','USER_SELECTED_ADJUSTED_REFERENCE',{'value':608,'units':'ohm/square','other_scenarios':[340,425]},'608 is model reference;340 room-temperature measurement;425 conditional RRR conversion, not an uncertainty interval.')
    setrow('sigma_n','DERIVED_JOINT_CONSTRAINT',{'value':sigma,'units':'S/m'},'1/(Rsheet*thickness); do not tune independently of D and N0.')
    setrow('N0','DERIVED_JOINT_CONSTRAINT',{'value':N0,'units':'J^-1 m^-3 per spin'},'sigma/(2 e^2 D); total-spin literature values are twice this convention.')
    setrow('gap_ratio','PRESERVED_MODEL_ASSUMPTION',{'value':1.7638769888620456,'units':'Delta0/(kB Tc)'},'Weak-coupling functional. Strong-coupling gap changes require consistent energy/current/frequency physics.')
    setrow('delta0','DERIVED_MODEL_VALUE',{'value':1.7638769888620456*kB*8.65,'units':'J'},'Not a measured gap for the selected film.')
    setrow('width','USER_SELECTED_REFERENCE',{'value':80,'units':'nm','other_samples':[60,100,120]},'Start from80nm; other widths identify different samples, not error bars.')
    setrow('thickness','REFERENCE_NOMINAL',{'value':7,'units':'nm'},'Nominal film thickness, not a spatial third dimension or zero-uncertainty measurement.')
    setrow('bath','REFERENCE_SCENARIO',{'value':.9,'units':'K'},'Bath of selected K20 comparison. Pilot theta=.12 is a different temperature.')
    setrow('lambda_L','DERIVED_OR_UNIDENTIFIED',None,'Do not inherit540nm independently of selected sigma/gap/current response. Recompute model response;64pH/square is a different experimental inductance estimate.')
    setrow('n_atom gphonon alpha2F material_preprocess rate_prefactor','ABSOLUTE_RATE_IDENTIFICATION_OPEN',None,'Conditional NbN shape retained; exact mode basis, unit conversion, normalization and finite support still need evidence. No guessed renormalization.')
    setrow('lambdaeph','SOURCE_SPECIFIC_NOT_UNIVERSAL_RANGE',None,'Synthetic lambda=.1, derived-file integral and DFT values describe different inputs; no universal NbN confidence bounds.')
    setrow('Egamma','DERIVED_FROM_SELECTED_WAVELENGTHS',{'wavelength_nm':[775,1550],'energy_eV':[1239.841984332/775,1239.841984332/1550]},'Nominal photon energy hc/lambda, not retained energy or detection efficiency.')
    setrow('handoff_time','USER_SELECTED_TIME_ORIGIN',{'simulation_origin':'t0=0','optical_delay_difference':None},'Time from transfer; wavelength-dependent omitted cascade delay remains unknown and separate.')
    setrow('spectral_injection','PRESERVED_MODEL_CLOSURE',{'rule':'D.30 phonon occupation proportional to alpha2F/gphonon; normalized to Ein'},'No simultaneous thermal-photon or delayed electronic source; spectral units and support must be admitted first.')
    setrow('fano event_ensemble jitter fit_energy_variances','DEFERRED_NOT_ACTIVATED',None,'First task is deterministic hotbelt formation and relative crossing latency. Published ensemble fits do not define new stochastic inputs.')
    setrow('Delta_state p_state n_state','DYNAMIC_STATE_CONSTRAINT',{'electron_occupancy':[0,1],'phonon_min':0},'Fields and populations evolve together; initial equilibrium must be recalculated for the new scenario, not copied from the pilot.')
    setrow('TE Tph_equiv','DERIVED_DIAGNOSTIC',None,'Equivalent temperature summarizes an energy moment and does not replace the distribution or uniquely determine forces.')
    setrow('delta_reg','PRESERVED_EFFECTIVE_PARAMETER_UNCALIBRATED',{'value':.1,'units':'Delta0','historical_sensitivity':[.05,.1,.2]},'Finite core completion, not mesh epsilon. Historical alternatives are sensitivity cases, not an admitted interval.')
    setrow('qbare Gamma K0 time_scale Qdelta DL Gface','DERIVED_JOINT_CONSTRAINT',None,'Compute from the selected functional, fields, populations and geometry; do not scan this as an independent material input.')
    setrow('D36','MANDATORY_CONTINUOUS_ADMISSIBILITY',{'condition':'minimum principal eigenvalue positive with sign resolved against numerical uncertainty'},'An indefinite principal symbol is a physical/mathematical failure, not fixed by smaller dt or clipping.')
    setrow('kwt_times','UNIDENTIFIED_EFFECTIVE_MOBILITY_WITH_REFERENCES',{'pair_definition':'tau_ee(Tc), tau_ep(Tc)','inherited_ps':[.5,2.47],'Allmaras_reference_ps':[5,24.7],'K20_reference_ps':[6,24.7]},'Distinct mobility scenarios; no KWT pair newly adopted. Factor10 is not a unit conversion. Do not relabel scattering times as energy-relaxation or BGK time.')
    setrow('taukin','NOT_IDENTIFIED',{'synthetic_ps':.7,'physical_range':None},'Energy-conserving redistribution closure. No reviewed measurement identifies its constant value; comparison to a microscopic ee operator remains necessary.')
    setrow('tauescape','REFERENCE_SCENARIOS_NOT_CONFIDENCE_INTERVAL',{'K20_reference_ps':20,'inherited_ps':15,'Allmaras_acoustic_ps':9.4,'Allmaras_other_fit_ps':80},'Phonon escape to substrate depends on film/interface. Do not interpret9.4–80ps as sample uncertainty.')
    setrow('mobility_floor','PRESERVED_CLOSURE',{'value':.9,'units':'K'},'Floor Tmob=max(TE,Tbath) uses the selected bath; not a second bath.')
    setrow('vector_potential','PRESERVED_APPROXIMATION',{'A':0},'No applied magnetic field or resolved self-field; restricted model assumption.')
    setrow('external_heating','INACTIVE_IN_PHOTON_SCENARIO',{'value':0},'Do not import synthetic stage2 heating into photon deposition.')
    setrow('heating_closure joule_allocation','PRESERVED_ENERGY_ACCOUNTING_REQUIRES_SPATIAL_CHECK',None,'Deposit normal Joule heat and QDelta once; preserve signed superconducting/reservoir work, not additional positive heat.')
    setrow('extra_phonon_physics','OMITTED_PHYSICS_LIMIT',None,'No lateral phonon transport or new phonon-phonon operator activated. Missing effects constrain interpretation, not just numerical precision.')
    setrow('face_rate legacy_thermal_catalog legacy_mesh','LEGACY_OR_SYNTHETIC_NOT_TRANSFERABLE',None,'Historical controls are not parameters of the new material/spatial kinetic discretization.')
    setrow('length_2D aspect_ratio resolved_length','PROPOSED_DOMAIN_DESIGN',{'initial_candidates_L_over_W':[4,6,8],'longer_window_candidates_L_over_W':[8,12]},'Conditional planning spans; full2D first until a1D interface is justified. Original1.5–6W was illustrative, not a hard cap.')
    setrow('observation_time','PROPOSED_OBSERVATION_WINDOW',{'initial_ps':50,'extended_ps':[100,200]},'No event by window end means censored/undetermined, not proved nondetection. Domain must be reassessed when window extends.')
    setrow('boundary_tolerance linear_diffusion_filter','ANALYTICAL_PREFILTER_NOT_ERROR_BOUND',{'illustrative_exterior_mass_fractions':[.01,.001]},'Free Gaussian instantaneous exterior mass. Actual admission is insensitivity of observables to artificial terminal placement.')
    setrow('lead_left lead_right interface_R boundary_distances geometry_extensions','CONDITIONAL_DOMAIN_OR_INTERFACE',None,'No independent lead-length interval before transverse dynamics and kinetic interface are verified. Record distances from deposition center separately from any threshold-defined region.')
    setrow('diffusion_scales','DERIVED_PLANNING_SCALES',{'ell0_nm':4.6984908958,'xi_c_nm':6.6446695476,'transverse_mode_ps':12.9691115062},'Normal linear reference only; not guaranteed hotbelt or homogenization times.')
    setrow('section_leads','DERIVED_GEOMETRY',{'value':80e-9*7e-9,'units':'m2'},'Same width and thickness if a uniform continuation is later justified.')
    setrow('insulating_sides reservoir_branch reservoir_population reservoir_load normal_trace interface_kinetic reservoir_work','REQUIRED_BOUNDARY_PHYSICS_NOT_YET_DYNAMICALLY_ADMITTED',None,'Use stationary branch consistent with circuit current and shared-energy kinetic populations; complete dynamic reservoir work and external normal trace before physical transients.')
    setrow('potential_reference','GAUGE_CHOICE',{'right_terminal_V':0},'Fixes the arbitrary potential reference, not an additional physical voltage source.')
    setrow('Rbias','PRESERVED_MEMORY_CIRCUIT_SCENARIO',{'value':1e4,'units':'ohm'},'Polarization resistance; no measured K20 tolerance found.')
    setrow('Lbias','PRESERVED_MEMORY_CIRCUIT_SCENARIO',{'value':1e-6,'units':'H'},'Bias-branch inductance; not detector kinetic inductance.')
    setrow('Rload','PRESERVED_MEMORY_CIRCUIT_SCENARIO',{'value':50,'units':'ohm'},'Load resistance at the defined readout port.')
    setrow('Ccouple','PRESERVED_MEMORY_CIRCUIT_SCENARIO',{'value':1e-10,'units':'F'},'100pF from memory CM; historical1pF key does not override the selected circuit.')
    setrow('Ltotal','PRESERVED_MEMORY_CIRCUIT_SCENARIO',{'value':1e-8,'units':'H','K20_additional_inductor_reference_H':96e-9},'10nH total memory reference and96nH added experimental inductor are different definitions. No silent substitution.')
    setrow('Lext','DERIVED_ENERGY_PARTITION',None,'Recompute total-minus-resolved inductance at a declared reference state when material/domain changes; do not reuse9.666757nH or subtract instantaneous hotspot inductance.')
    setrow('Ib Is Vbias','PROPOSED_REFERENCE_BIAS_POINTS',{'initial_current_A':[15.5e-6,21.5e-6],'bias_voltage_V':[.155,.215]},'Published80nm comparison points; Vbias=Rbias Iinitial. Ib and Is subsequently evolve by CM and are not clamped. Require a stable stationary branch; not a certified continuous current interval.')
    setrow('vc','INITIAL_EQUILIBRIUM_CONDITION',{'initial_V':0},'Stationary superconducting DC state; capacitor voltage evolves after the event.')
    setrow('Vdev Vout latency observables','SEPARATE_OBSERVABLE_DEFINITIONS',None,'Device power-port voltage, load voltage and condensate events are distinct. Relative deterministic crossing latency is not an IRF peak difference.')
    setrow('threshold readout_windows recovery electronics_transfer noise optical_chain tapers experimental_readout observation_uncertainty','EXPERIMENTAL_MAPPING_INCOMPLETE',None,'Memory observable conventions remain identifiable, but no complete calibrated K20 transfer/threshold/noise/position mapping is supplied. Do not invent a physical range from software defaults.')
    setrow('amplitude_box gamma_box','AVAILABLE_TABLE_COVERAGE_ONLY',{'positive_amplitude_Delta0':[.08,1.5],'Gamma_Delta0':[0,1.2]},'Not a Cartesian physical confidence domain. The separate exact normal point does not validate the0-to-.08 gap.')
    setrow('eta','RETAINED_VERIFIED_NUMERICAL_CONTROL',{'value':1e-8,'units':'Delta0'},'Causal regulator, not material broadening; scope of convergence must follow the used state.')
    setrow('xmax electron_nodes phonon_grid OmegaD count_partition reaction_quadrature face_quadrature','REQUIRES_MATERIAL_SUPPORT_REDESIGN',{'electronic_candidate_nodes':630,'phonon_candidate_nodes':1025,'conditional_NbN_positive_support_max_meV':67.784,'sampled_electronic_max_meV':15.832,'synthetic_phonon_max_meV':5.259},'See material_kinetics.md section7. xmax=12 is a count coordinate, not exactly Emax/Delta0. Synthetic support and resolution are not NbN rates; redesign extent, native energy, tails and event support before setting counts.')
    setrow('matsubara newton','ALGORITHMIC_CONTROL_NOT_PHYSICAL_RANGE',None,'Retain causal solve and energy accuracy; no arbitrary relaxation of root tolerance to match experimental uncertainty.')
    setrow('predictor','LOCAL_ACCELERATOR_ONLY',{'amplitude_Delta0':[.985,1.005],'Gamma_Delta0':[.003,.009]},'Existing seed box does not cover hotbelt. Design direct evaluation or explicitly extend/reverify; no silent extrapolation.')
    setrow('sem_degree','PROPOSED_NUMERICAL_START',{'value':4},'Maintain spatial representation used in stage3; degree does not by itself resolve a core.')
    setrow('spatial_mesh','PROPOSED_RESOLUTION',{'smooth_element_nm':[10,20],'core_comparison_element_nm':[5,10],'conditional_compact_deposition_local_element_nm':[2.1,2.9]},'Element length, not nodal separation. Largest degree4GLL gap=.3273h. Compact s1.4–1.9nm requires local refinement; h<=s/(2*.3273) is only a starting resolution heuristic, not an accuracy certificate. Graded2D mesh should avoid uniform cost.')
    setrow('fd_controls static_solver','DIAGNOSTIC_SPECIFIC_NOT_PHYSICAL_MARGIN',None,'Historical derivative increments and minimizer tolerances keep their original scope; stationary force/current must be controlled for the new initialization.')
    setrow('time_method dt','UNADMITTED_SPATIAL_INTEGRATOR_DESIGN',{'proposed_initial_step_fs':[1,10]},'Not a stable-step guarantee; coupled spectrum, mesh and collision rates determine admissible steps. No physical method chosen merely by literature.')
    setrow('checkpoint','PROPOSED_OUTPUT_CADENCE',{'early_field_spacing_ps':[.05,.1]},'Saved field cadence is distinct from solver time step and event location accuracy.')
    setrow('numerical_budgets','PROPOSED_OBSERVABLE_BUDGET',{'relative_crossing_latency_difference_ps':.1,'exploratory_ps':.2},'Absolute numerical sensitivity target for the difference, not experimental CI or dt. Do not replace positivity/support/balance with this target.')
    setrow('stopping','STRUCTURAL_STOPS_PRESERVED',None,'Stop and report unsupported populations, lost causality, nonpositive D36 or invalid branch; never clip to force acceptance.')
    setrow('Ic_calibration','NO_TARGET_REFIT_BY_DEFAULT',None,'Do not retune D to an inherited critical current while claiming independent experimental diffusion; test selected points against actual branch.')
    setrow('material_provenance run_identity execution_resources','TRACEABILITY_AND_COST_RULE',{'agent_expected_max_seconds':300,'uncertain_timeout_seconds':240},'Immutable evidence and fresh directories. Long commands only in Geminga notebook plus chat; no heavy run executed in this research.')
    setrow('sample_identity substrate phonon_gamma_legacy','SOURCE_SPECIFIC_REFERENCE',None,'Selected80nm K20 NbN/SiO2/Si case; legacy phonon heat-capacity gamma is not Usadel pair-breaking Gamma or an admitted full-spectrum kernel.')
    # User decision is read explicitly: absence never means approval.
    decisions=json.loads((OUT/'user_decisions.json').read_text(encoding='utf8'))
    transfer=decisions['source_based_transfer_proposal']
    setrow('retention','SOURCE_REFERENCE_PROPOSED',transfer['retention'],
        'K20 fitted reference fraction, not a standard universal measurement. Translation to D.30 must not include a second subtraction of post-transfer losses.')
    setrow('Ein','SOURCE_REFERENCE_PROPOSED',transfer['Ein_eV'],
        'Conditional retained energies in eV using hc/lambda and the K20 fraction0.667; not optical absorption/detection efficiency.')
    setrow('photon_sigma',transfer['gaussian_sigma_nm']['status'],transfer['gaussian_sigma_nm'],
        'Gaussian standard deviation per coordinate. A20 total-energy R90 bridge, not a measured Korzh phonon width or validated interval. Depth averaging and electron residual remain explicit.')
    setrow('photon_position','CONTROLLED_REFERENCE_PROPOSED',transfer['position'],
        'Central controlled deposition tests formation, not an imposed hotbelt or experimental absorption-position ensemble.')
    entries=[]
    for p in inventory['parameters']:
        doc,default=domains[p['family']]
        status,value,note=specific.get(p['id'],(default,None,'No independent interval identified; consult the family review and original source-specific constraints.'))
        entries.append(dict(id=p['id'],family=p['family'],symbol=p['symbol'],physical_name=p['name'],units=p['units'],
            status=status,value_or_proposal=value,rationale=note,review=doc,physical_confidence_interval=None,
            original_source_refs=p['source_refs']))
    assert len(entries)==127 and len({x['id'] for x in entries})==127
    assert not set(specific)-{x['id'] for x in entries}
    missing=[p['id'] for p in inventory['parameters'] if p['id'] not in specific]
    result=dict(schema='pysnspd.stage3_5.range_ledger.r1',status='RESEARCH_WITH_PROVISIONAL_SCENARIOS_AND_EXPLICIT_UNKNOWN_RANGES',
      base_inventory_sha256=hashlib.sha256(BASE.read_bytes()).hexdigest(),user_decisions='user_decisions.json',
      parameters=entries,count=len(entries),families=len(domains),independent_physical_confidence_box_admitted=False,
      code_parameters_changed=False,new_physical_transients=0,unresolved_family_default_ids=missing)
    (OUT/'range_ledger.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print(json.dumps(dict(count=len(entries),specific_decisions=len(specific),family_only=missing,scenario_sigma=sigma,N0_per_spin=N0)))

if __name__=='__main__':main()
