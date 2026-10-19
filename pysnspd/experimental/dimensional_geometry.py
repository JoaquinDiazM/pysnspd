"""Deterministic convex-domain repair of a malformed boundary Voronoi cell.

Keep sites/triangles/edges of the inherited Delaunay mesh. Reconstruct exact
clipped Voronoi volumes and common faces by neighboring-site half planes.
Never renormalize all volumes to hide a missing cell or search mesh seeds.
"""
import json,hashlib
from pathlib import Path
import numpy as np
from .effective_bcs_preparation import load_dual_geometry


def clipped_voronoi(coordinates_m,edges,length_m,width_m):
    xy=np.asarray(coordinates_m,float)*1e9
    L,W=length_m*1e9,width_m*1e9
    adjacency=[[] for _ in xy]
    for i,j in edges:adjacency[i].append(j);adjacency[j].append(i)
    polygons=[];areas=[]
    def clip(poly,normal,offset):
        out=[]
        for start,end in zip(poly,np.roll(poly,-1,axis=0)):
            a=float(start@normal-offset);b=float(end@normal-offset)
            ina=a<=0.;inb=b<=0.
            if ina:out.append(start)
            if ina!=inb:out.append(start+(end-start)*(a/(a-b)))
        return np.asarray(out).reshape(-1,2)
    for i,site in enumerate(xy):
        poly=np.array([[0.,-W/2],[L,-W/2],[L,W/2],[0.,W/2]])
        for j in adjacency[i]:
            normal=xy[j]-site
            offset=.5*np.dot(xy[j]+site,normal)
            poly=clip(poly,normal,offset)
            if len(poly)<3:raise ArithmeticError("Empty clipped Voronoi cell")
        # Translate before the shoelace sum to avoid subtracting large x^2 terms.
        local=poly-site
        area=.5*abs(np.dot(local[:,0],np.roll(local[:,1],-1))-
                      np.dot(local[:,1],np.roll(local[:,0],-1)))
        polygons.append(poly);areas.append(area*1e-18)
    faces=[]
    for i,j in edges:
        normal=xy[j]-xy[i];offset=.5*np.dot(xy[j]+xy[i],normal)
        poly=polygons[i]
        tolerance=2048*np.finfo(float).eps*max(1.,abs(offset),np.linalg.norm(normal)*L)
        on=poly[abs(poly@normal-offset)<=tolerance]
        if len(on)<2:faces.append(0.);continue
        tangent=np.array([-normal[1],normal[0]])/np.linalg.norm(normal)
        s=on@tangent
        faces.append(float(np.ptp(s))*1e-9)
    areas=np.asarray(areas);faces=np.asarray(faces)
    if np.any(areas<=0) or not np.isclose(areas.sum(),length_m*width_m,rtol=2e-11,atol=0):
        raise ArithmeticError("Clipped Voronoi cells do not exactly cover strip")
    return areas,faces


def verified_geometry(args,output,ell0_m,legacy):
    output=Path(output);target=output/"meshes/v12_boundary_repair/mesh.npz"
    receipt=target.with_suffix(".json")
    if target.exists():
        info=json.loads(receipt.read_text())
        digest=hashlib.sha256(target.read_bytes()).hexdigest()
        if digest!=info["mesh_sha256"]:raise ValueError("Verified v1.2 mesh changed")
        result=load_dual_geometry(target,expected_sha256=digest)
        if not np.isclose(result.ell0_m,ell0_m,rtol=1e-13,atol=0):
            raise ValueError("Verified mesh normalization changed")
        result.receipt.update(info);return result
    try:return legacy(args,output,ell0_m)
    except ValueError as exc:
        if str(exc)!="Dual areas do not cover the physical strip":raise
    raw=output/"meshes/reference_geometry"/(
        "l"+format(args.length_nm,"g")+"_w80_h"+format(args.spacing_nm,"g"))/"mesh.npz"
    with np.load(raw,allow_pickle=False) as z:
        data={key:z[key].copy() for key in z.files}
    original=hashlib.sha256(raw.read_bytes()).hexdigest()
    xy=data["coordinates_m"];edges=data["edges"]
    areas,faces=clipped_voronoi(xy,edges,args.length_nm*1e-9,80e-9)
    lengths=np.linalg.norm(xy[edges[:,1]]-xy[edges[:,0]],axis=1)
    data.update(dual_area_m2=areas,area_weights=areas/ell0_m**2,
        dual_face_length_m=faces,conductance=faces/lengths,
        coordinates_bar=(xy-data.get("origin_m",np.array([args.length_nm*.5e-9,0.])))/ell0_m,
        ell0_m=np.array(ell0_m))
    target.parent.mkdir(parents=True,exist_ok=True)
    with target.open("wb") as f:np.savez_compressed(f,**data)
    digest=hashlib.sha256(target.read_bytes()).hexdigest()
    info=dict(mesh_sha256=digest,original_mesh_sha256=original,
        original_mesh_path=str(raw),boundary_repair="Clipped convex Voronoi half planes",
        triangulation_unchanged=True,sites_unchanged=True,no_global_area_renormalization=True,
        repaired_area_relative_error=float(areas.sum()/(args.length_nm*1e-9*80e-9)-1))
    receipt.write_text(json.dumps(info,indent=2)+"\n")
    result=load_dual_geometry(target,expected_sha256=digest);result.receipt.update(info)
    return result
