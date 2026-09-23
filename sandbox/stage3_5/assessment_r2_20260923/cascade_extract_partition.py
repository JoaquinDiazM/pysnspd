"""Digitize A20 Figs2.6(a)/2.8(a), not a cascade simulation.

Use embedded original RGB bytes, visible solid/dashed/dotted curves only.
At a finite radius, positivity and no escape imply
  E_e(R)/Egamma <= E_e(infinity)/Egamma <= 1-E_ph(R)/Egamma.
No terminal sample is silently normalized to one.
"""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path
from statistics import median
import pdfplumber
from PIL import Image

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT/"tmp/pdfs/modelo_v0_3/sources/Allmaras_Thesis_Final.pdf"
OUT = ROOT/"docs/implementation/stage3_5/assessment_r2_20260923/cascade/cascade_partition.json"
EXPECTED = "2554e9a2b047f93f7c7f50a2897887afaf522eb530e057307c5a31ef4e1e0e96"
A6 = {"left":252., "right":1583., "top":188., "bottom":1141., "xmax":10.,
      "legend":(962,294,1552,1038)}
A8 = {"left":252., "right":1582., "top":139., "bottom":1083., "xmax":2.,
      "legend":(1080,174,1553,645)}
C6 = {"t5e-5":(204,121,15), "t1e-4":(204,0,0)}
C8 = {2:(0,31,170), 3:(0,163,20), 4:(166,163,20), 5:(204,0,0)}


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def embedded(page, idx):
    obj=page.images[idx]
    data=obj["stream"].get_data()
    im=Image.frombytes("RGB",tuple(obj["srcsize"]),data)
    return im, {"name":obj["name"],"size_px":list(obj["srcsize"]),
                "RGB_sha256":hashlib.sha256(data).hexdigest(),
                "vector_curves":len(page.curves),"vector_lines":len(page.lines)}


def frac(y, axes):
    return (axes["bottom"]-y)/(axes["bottom"]-axes["top"])


def visible(x,y,axes):
    a,b,c,d=axes["legend"]
    return not (a<=x<=c and b<=y<=d)


def point_cloud(im,axes,color,tolerance=12):
    px=im.load()
    result={}
    for x in range(int(axes["left"])+3,int(axes["right"])-3):
        ys=[y for y in range(int(axes["top"])+3,int(axes["bottom"])-3)
            if visible(x,y,axes) and
            max(abs(v-c) for v,c in zip(px[x,y],color))<=tolerance]
        if ys:
            result[x]=ys
    return result


def selected(cloud,axes,low,high):
    return {x:[y for y in ys if low<frac(y,axes)<high] for x,ys in cloud.items()
            if any(low<frac(y,axes)<high for y in ys)}


def sample(cloud,axes,xvalue):
    x0=axes["left"]+xvalue/axes["xmax"]*(axes["right"]-axes["left"])
    samples=[(x,median(ys),min(ys),max(ys)) for x,ys in cloud.items() if abs(x-x0)<=12]
    if len(samples)<4:
        raise ValueError(f"No visible local trace at {xvalue}")
    xm=sum(v[0] for v in samples)/len(samples)
    ym=sum(v[1] for v in samples)/len(samples)
    slope=sum((x-xm)*(y-ym) for x,y,*_ in samples)/sum((x-xm)**2 for x,*_ in samples)
    y0=ym+slope*(x0-xm)
    residual=max(abs(y-(y0+slope*(x-x0))) for x,y,*_ in samples)
    halfstroke=max((hi-lo)/2 for _,_,lo,hi in samples)
    padding=residual+halfstroke+2.5+abs(slope)*2
    values=[(bottom-y)/(bottom-top)
            for bottom in (axes["bottom"]-2,axes["bottom"]+2)
            for top in (axes["top"]-2,axes["top"]+2)
            for y in (y0-padding,y0+padding)]
    return {"fraction_of_photon":frac(y0,axes),
            "reading_envelope":[min(values),max(values)],
            "sample_center_px":[x0,y0],"local_slope":slope,
            "stroke_halfwidth_px":halfstroke,"fit_residual_px":residual,
            "sample_count":len(samples),
            "sampled_centers_px":[[x,y] for x,y,*_ in samples]}


def inverse_visible(cloud,axes,target):
    samples=sorted((x,frac(median(ys),axes)) for x,ys in cloud.items())
    crossings=[]
    for (x0,f0),(x1,f1) in zip(samples,samples[1:]):
        if 0<x1-x0<=28 and min(f0,f1)<=target<=max(f0,f1) and f0!=f1:
            crossings.append(x0+(target-f0)*(x1-x0)/(f1-f0))
    if not crossings:
        raise ValueError(f"No visible inverse crossing for fraction {target}")
    # This inversion is performed only on the positive monotone CDF segments.
    # Small pixel-level backsteps may give several crossings; keep their extent.
    values=[(x-axes["left"])/(axes["right"]-axes["left"])*axes["xmax"] for x in crossings]
    return median(values),[min(values),max(values)]


def phonon_radius(cloud,axes,q,global_ph_bounds):
    # Radii bracket all target normalizations and a conservative pointwise
    # vertical reading budget: 6px / plot height plus +/-2px horizontal axes.
    vertical=6/(axes["bottom"]-axes["top"])
    lo=q*global_ph_bounds[0]-vertical
    hi=q*global_ph_bounds[1]+vertical
    low,_=inverse_visible(cloud,axes,lo)
    high,_=inverse_visible(cloud,axes,hi)
    hpad=4*axes["xmax"]/(axes["right"]-axes["left"])
    return {"quantile":q,"radius_envelope_nm":[low-hpad,high+hpad],
            "global_phonon_fraction_used":global_ph_bounds,
            "CDF_target_fraction_of_photon_envelope":[q*global_ph_bounds[0],q*global_ph_bounds[1]],
            "vertical_reading_budget_fraction":vertical,
            "normalization":"own phonon energy, bounded from positivity and no escape"}


def main():
    if digest(SOURCE)!=EXPECTED:
        raise ValueError("PDF changed: re-inspect image calibration")
    with pdfplumber.open(SOURCE) as pdf:
        im6,m6=embedded(pdf.pages[33],0)
        im8,m8=embedded(pdf.pages[34],2)
    rows=[]
    cache6={}
    for key,color in C6.items():
        cloud=point_cloud(im6,A6,color)
        e=selected(cloud,A6,0.,.14)
        ph=selected(cloud,A6,.3,1.)
        cache6[key]=(e,ph)
        e95,ph95=sample(e,A6,9.5),sample(ph,A6,9.5)
        elo,ehi=e95["reading_envelope"]
        plo,phi=ph95["reading_envelope"]
        eg=[elo,1-plo]
        pg=[plo,1-elo]
        if not 0<eg[0]<eg[1]<1 or not 0<pg[0]<pg[1]<1:
            raise ValueError("Invalid finite-radius global bound")
        r50=phonon_radius(ph,A6,.5,pg)
        r90=phonon_radius(ph,A6,.9,pg)
        ratio=[r90["radius_envelope_nm"][0]/r50["radius_envelope_nm"][1],
               r90["radius_envelope_nm"][1]/r50["radius_envelope_nm"][0]]
        nominal_tau=5e-5 if key=="t5e-5" else 1e-4
        rows.append({"key":key,"time_over_tau0":nominal_tau,
                     "time_ps":nominal_tau*1870,
                     "outer_radius_nm":9.5,
                     "electron_inside_outer_radius":e95,
                     "phonon_inside_outer_radius":ph95,
                     "electron_global_fraction_bound":eg,
                     "phonon_global_fraction_bound":pg,
                     "outside_outer_radius_energy_fraction_bound":[0.,1-elo-plo],
                     "phonon_R50":r50,"phonon_R90":r90,
                     "phonon_R90_over_R50_envelope":ratio,
                     "single2DGaussian_ratio_compatible":
                         ratio[0]<=math.sqrt(math.log(10)/math.log(2))<=ratio[1],
                     "at_radius5":{"electron":sample(e,A6,5.),
                                   "phonon":sample(ph,A6,5.)}})
    fig8=[]
    # t<1.2 on this axis avoids the legend completely for all three channels.
    for radius,color in C8.items():
        cloud=point_cloud(im8,A8,color)
        for xvalue in (.3,.5,1.):
            center=A8["left"]+xvalue/2*(A8["right"]-A8["left"])
            nearby=sorted(set(y for x,ys in cloud.items() if abs(x-center)<=12 for y in ys))
            groups=[]
            for y in nearby:
                if not groups or y-groups[-1][-1]>7:
                    groups.append([y])
                else:
                    groups[-1].append(y)
            if len(groups)!=3:
                raise ValueError(f"Ambiguous total/phonon/electron branches: r={radius},t={xvalue},groups={groups}")
            measurements={}
            for label,g in zip(("total","phonon","electron"),groups):
                branch={x:[y for y in ys if min(g)<=y<=max(g)] for x,ys in cloud.items()
                        if abs(x-center)<=12 and any(min(g)<=y<=max(g) for y in ys)}
                measurements[label]=sample(branch,A8,xvalue)
            residual=measurements["total"]["fraction_of_photon"]-sum(
                measurements[k]["fraction_of_photon"] for k in ("electron","phonon"))
            fig8.append({"radius_nm":radius,"time_axis_value":xvalue,
                         "time_over_tau0":xvalue*1e-4,"time_ps":xvalue*1e-4*1870,
                         **measurements,"total_minus_sum_reading_residual":residual})
    checks=[]
    for row in rows:
        other=next(r for r in fig8 if r["radius_nm"]==5 and
                   abs(r["time_ps"]-row["time_ps"])<1e-12)
        for label in ("electron","phonon"):
            a=row["at_radius5"][label]
            b=other[label]
            bounds=[max(a["reading_envelope"][0],b["reading_envelope"][0]),
                    min(a["reading_envelope"][1],b["reading_envelope"][1])]
            checks.append({"time_ps":row["time_ps"],"radius_nm":5,"channel":label,
                           "difference_fig6_minus_fig8":a["fraction_of_photon"]-b["fraction_of_photon"],
                           "envelopes_overlap":bounds[0]<=bounds[1]})
    result={
       "schema":"pysnspd.stage3_5.cascade_partition_digitization.v1",
       "status":"HISTORICAL_A20_REFERENCE_ONLY_WIDTH_OPEN","new_physics_simulations":0,
       "source":{"path":str(SOURCE.relative_to(ROOT)).replace("\\","/"),
                 "sha256":EXPECTED,"url":"https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf",
                 "fig2_6a":{"PDF_page":34,"printed_page":22,**m6},
                 "fig2_8a":{"PDF_page":35,"printed_page":23,**m8}},
       "script_sha256":digest(Path(__file__)),
       "scenario":{"photon_energy_eV":1,"Tb_K":4.325,"D_cm2_s":.5,
                   "geometry":"cylindrical","initial_electron_pairs":2,
                   "imposed_initial_radius_nm":1,"tau0_kinetic_ps":1870,
                   "ee_phonon_diffusion_escape":"omitted in source radial comparison",
                   "normalization":"all raw plotted energies divided by photon energy"},
       "axes":{"fig2_6a":A6,"fig2_8a":A8},
       "reading_policy":{"RGB_Linf_tolerance":12,"axis_position_uncertainty_px":2,
                         "figure6_inverse_vertical_budget_px":6,
                         "figure6_max_dashed_interpolation_gap_px":28,
                         "legend_regions_excluded":True,
                         "no_endpoint_extrapolation":True,
                         "envelopes":"deterministic reading budgets, not statistical or physical confidence"},
       "finite_radius_and_phonon_percentiles":rows,
       "figure2_8_visible_samples":fig8,
       "cross_figure_checks":checks,
       "gaussian2D_R90_over_R50":math.sqrt(math.log(10)/math.log(2)),
       "restrictions":["Finite-radius global bounds assume excess nonnegative and exactly no source-model escape",
                       "No Gaussian width adopted","No rescaling to Korzh colors or bath",
                       "No density peak, complete second moment or spectrum extracted",
                       "No full3D surfaces or absorption-depth distribution identified"]}
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"output":str(OUT.relative_to(ROOT)),
      "rows":[{k:r[k] for k in ("time_ps","electron_global_fraction_bound",
         "phonon_global_fraction_bound","phonon_R50","phonon_R90",
         "phonon_R90_over_R50_envelope","single2DGaussian_ratio_compatible")} for r in rows],
      "cross_figure_checks":checks,
      "max_total_identity_reading_residual":max(abs(r["total_minus_sum_reading_residual"]) for r in fig8)},
      indent=2))


if __name__=="__main__":
    main()
