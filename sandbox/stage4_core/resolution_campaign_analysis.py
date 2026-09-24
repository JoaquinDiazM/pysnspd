"""Compare registered moment variations without inventing a DOS tolerance."""
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator

OBSERVABLES = {
    'charge_current':('current','edge'),
    'energy_weighted_energy_flux':('energy_flux','edge'),
    'gap_force_density':('force','node_integrated'),
    'amplitude_force_density':('amplitude_force_density','node'),
    'phase_torque_density':('phase_torque_density','node'),
}


def rectangular_interpolate(coordinates, values, targets):
    coordinates, targets = np.asarray(coordinates), np.asarray(targets)
    xs, ys = np.unique(coordinates[:, 0]), np.unique(coordinates[:, 1])
    if len(xs)*len(ys) != len(coordinates):
        raise ValueError('Registered interpolation requires a complete rectangular grid')
    grid = np.empty((len(xs), len(ys)), dtype=np.asarray(values).dtype)
    grid[np.searchsorted(xs,coordinates[:, 0]), np.searchsorted(ys,coordinates[:, 1])] = values
    return RegularGridInterpolator((xs,ys), grid, bounds_error=True)(targets)


def transfer_flux(source, target, current):
    """Interpolate J/dual-width at edge centers, then integrate target faces.

    This avoids comparing mesh-dependent raw link currents or mixing horizontal
    and vertical components. It is an observed interpolation comparison, not a
    certified continuum error bound.
    """
    sx, tx = source['coordinates_bar'], target['coordinates_bar']
    se, te = source['edges'], target['edges']
    sv, tv = sx[se[:, 1]]-sx[se[:, 0]], tx[te[:, 1]]-tx[te[:, 0]]
    sl, tl = np.linalg.norm(sv,axis=1), np.linalg.norm(tv,axis=1)
    sw, tw = source['conductance']*sl, target['conductance']*tl
    if np.any(sw<=0) or np.any(tw<=0):
        raise ValueError('Registered rectangular comparison needs positive dual widths')
    sm, tm = (sx[se[:, 0]]+sx[se[:, 1]])/2, (tx[te[:, 0]]+tx[te[:, 1]])/2
    result = np.empty(len(te),dtype=np.asarray(current).dtype)
    for axis in (0,1):
        select = abs(sv[:,axis])>abs(sv[:,1-axis])
        wanted = abs(tv[:,axis])>abs(tv[:,1-axis])
        density = current[select]/sw[select]*np.sign(sv[select,axis])
        sampled = rectangular_interpolate(sm[select],density,tm[wanted])
        result[wanted] = sampled*tw[wanted]*np.sign(tv[wanted,axis])
    return result


def observable(fields, probe, state, name):
    suffix, kind = OBSERVABLES[name]
    value = fields[probe+'_'+state+'_'+suffix]
    return value/fields['area_weights'] if kind=='node_integrated' else value


def map_observable(source, target, value, name):
    if np.array_equal(source['coordinates_bar'],target['coordinates_bar']):
        return value
    kind = OBSERVABLES[name][1]
    return (transfer_flux(source,target,value) if kind=='edge' else
            rectangular_interpolate(source['coordinates_bar'],value,target['coordinates_bar']))


def geometry_weights(fields, name, radius):
    if OBSERVABLES[name][1]=='edge':
        xy,edges,c = fields['coordinates_bar'],fields['edges'],fields['conductance']
        midpoint = (xy[edges[:,0]]+xy[edges[:,1]])/2
        mask=(c>0)&(np.linalg.norm(midpoint,axis=1)<=radius)
        weights=np.zeros(len(c));weights[mask]=1/c[mask]
    else:
        weights=fields['area_weights'].copy()
        weights[np.linalg.norm(fields['coordinates_bar'],axis=1)>radius]=0
        weights[fields['boundary_nodes']]=0
    return weights


def metrics(fields, probe, name, radius):
    weights=geometry_weights(fields,name,radius)
    norm=lambda value:float(np.sqrt(np.dot(weights,abs(value)**2)))
    full=observable(fields,probe,'full',name)
    projected=observable(fields,probe,'projected',name)
    zero=observable(fields,probe,'zero_hT',name)
    D,C=norm(projected-full),norm(zero-full)
    pf,ff=norm(projected),norm(full)
    alignment=float(np.real(np.vdot(full*weights,projected))/(pf*ff)) if pf*ff>1e-24 else None
    return dict(projection_difference_D=D,omitted_response_C=C,full_norm=ff,projected_norm=pf,
        D_over_C=(D/C if C>1e-12 else None),weighted_alignment=alignment,
        tiny_denominator_policy='Report absolute norms when a dimensionless norm is <=1e-12; this is presentation, not an acceptance tolerance')


def analyze(plan, output, projection):
    fields={}
    rows={row['id']:row for row in projection['records']}
    for key,row in rows.items():
        with np.load(Path(output)/'projection'/row['fields_path']) as value:
            fields[key]={name:value[name].copy() for name in value.files}
    radius=plan['kinetic']['probe_radius_ell0']
    moments={key:{probe:{name:metrics(value,probe,name,radius) for name in OBSERVABLES}
                  for probe in plan['kinetic']['probes']} for key,value in fields.items()}
    comparisons={}
    def compare(label,reference,alternative):
        target,source=fields[reference],fields[alternative]
        result={}
        for probe in plan['kinetic']['probes']:
            result[probe]={}
            for name in OBSERVABLES:
                weights=geometry_weights(target,name,radius)
                result[probe][name]={}
                for state in ('full','projected','zero_hT'):
                    a=observable(target,probe,state,name)
                    b=map_observable(source,target,observable(source,probe,state,name),name)
                    result[probe][name][state]=float(np.sqrt(np.dot(weights,abs(a-b)**2)))
        comparisons[label]=dict(reference=reference,alternative=alternative,difference_norms=result,
            interpretation='Observed field difference in the same weighted norm; not a rigorous bound')
    for group in plan['groups']:
        if group['grid']=='fine50':
            compare('energy_'+group['id'],group['id']+'_fine50',group['id']+'_base31')
    ids={(g['case_id'],g['eta_relative']):g['id'] for g in plan['groups']}
    reference=ids['radial_65_N256',.01]+'_fine50'
    compare('contour_eta',reference,ids['radial_65_N256',.02]+'_fine50')
    compare('spatial_mesh',reference,ids['radial_129_N256',.01]+'_fine50')
    budgets={}
    relevant=['energy_'+ids['radial_65_N256',.01],'contour_eta','spatial_mesh']
    for probe in plan['kinetic']['probes']:
        budgets[probe]={}
        for name in OBSERVABLES:
            terms=[comparisons[key]['difference_norms'][probe][name] for key in relevant]
            U_D=sum(term['full']+term['projected'] for term in terms)
            U_C=sum(term['full']+term['zero_hT'] for term in terms)
            value=moments[reference][probe][name];D,C=value['projection_difference_D'],value['omitted_response_C']
            if D<=U_D:
                assessment='Projection discrepancy unresolved against observed numerical variation'
            elif D>=C+U_D+U_C:
                assessment='Projection does not improve this moment beyond observed variation'
            elif D+U_D<C-U_C:
                assessment='Projection improves this moment beyond observed variation; not dynamic admission'
            else:
                assessment='Projection discrepancy resolved; improvement over omitted response remains uncertain'
            budgets[probe][name]=dict(**value,observed_variation_U_D=U_D,observed_variation_U_C=U_C,
                components=relevant,assessment=assessment,
                alignment_warning=(value['weighted_alignment'] is not None and value['weighted_alignment']<0))
    return dict(schema='pysnspd.stage4.moment_resolution_analysis.v1',moments=moments,
        comparisons=comparisons,reference_with_all_three_controls=reference,diagnostic_variation_budget=budgets,
        chi_quadrature={key:dict(numerical=row['chi_quadrature'],exact=row['chi_integral_exact_window']) for key,row in rows.items()},
        conclusion_policy='Interpret moment differences relative to observed grid, contour and mesh changes; no fixed DOS tolerance and no automatic stage closure.',
        mesh_transfer='Interpolate nodal force DENSITIES and edge flux per dual width; then evaluate on the65x65 geometry.',
        uncertainty_scope='U is a sum of observed full/projected or full/zero changes, not a certified error bound. Asymmetric base31 has no independent fine-grid control.',
        physical_time_steps=0,stage4_complete=False,production_changed=False)
