"""Parallel numerical acceleration of the native BCS event algebra.

No physical closure, rates or support policy differs from the NumPy oracle.
Numba cache is disabled: no builds are written into the source repository.
"""
import math
import numpy as np
from numba import njit, prange, set_num_threads, get_num_threads


@njit(cache=False,inline='always')
def _difference(a,b):
    if a>b:return math.exp(a)*(-math.expm1(b-a))
    if b>a:return -math.exp(b)*(-math.expm1(a-b))
    return 0.


@njit(cache=False,inline='always')
def _interp_log(logrow,lo,hi,b):
    out=0.
    if b<1:out+=(1-b)*logrow[lo]
    if b>0:out+=b*logrow[hi]
    return out


@njit(cache=False,parallel=True)
def acoustic_collisions(e,amplitude,p,n,lp,lh,ln,l1n,i,j,rec,pairweights,
                         phenergy,phcapacity,w,coefficient,cutlow,cuthigh,prefactor,equilibrium_mask):
    nn,ne=e.shape;npop=n.shape[1];nev=len(i)
    rate=np.zeros((nn,nev));omega=np.zeros((nn,nev));beta=np.zeros((nn,nev))
    lower=np.zeros((nn,nev),np.int64);upper=np.ones((nn,nev),np.int64)
    dp=np.zeros_like(p);dn=np.zeros_like(n);invalid=np.zeros(nn,np.int64)
    for node in prange(nn):
        if equilibrium_mask[node]:continue
        a2=amplitude[node]**2
        for ev in range(nev):
            ii,jj=i[ev],j[ev]
            ei,ej=e[node,ii],e[node,jj]
            om=ei+ej if rec[ev] else ej-ei
            omega[node,ev]=om
            if om<cutlow or om>cuthigh:continue
            up=max(1,np.searchsorted(phenergy,om))
            if up>=npop or om<phenergy[0]:
                invalid[node]=1
                continue
            lo=up-1;b=(om-phenergy[lo])/(phenergy[up]-phenergy[lo])
            lower[node,ev]=lo;upper[node,ev]=up;beta[node,ev]=b
            coh=1+a2/(ei*ej) if rec[ev] else 1-a2/(ei*ej)
            c=prefactor*pairweights[ev]*coefficient*om*om*coh
            emit=_interp_log(l1n[node],lo,up,b)
            absorb=_interp_log(ln[node],lo,up,b)
            if rec[ev]:
                f=lp[node,ii]+lp[node,jj]+emit
                r=lh[node,ii]+lh[node,jj]+absorb
            else:
                f=lh[node,ii]+lp[node,jj]+emit
                r=lp[node,ii]+lh[node,jj]+absorb
            value=c*_difference(f,r);rate[node,ev]=value
            dp[node,ii]+= -value if rec[ev] else value
            dp[node,jj]-=value
            dn[node,lo]+=(1-b)*value
            dn[node,up]+=b*value
        for k in range(ne):dp[node,k]/=4*w[k]
        for k in range(npop):dn[node,k]/=phcapacity[k]
    if np.any(invalid):raise ValueError('Coupled event lacks phonon bracket; extend physical quadrature')
    return dp,dn,rate,omega,lower,upper,beta


@njit(cache=False)
def transport(e,amplitude,lp,lh,xi,w,edges,volumes,face_rates,x,gw,N0,Eu):
    nn,ne=e.shape;nf=len(edges);nq=len(x)
    flux=np.zeros((nf,nq));q=np.zeros((nf,nq))
    lower=np.zeros((2,nf,nq),np.int64);upper=np.zeros((2,nf,nq),np.int64)
    beta=np.zeros((2,nf,nq));lpi=np.zeros((2,nf,nq));lhi=np.zeros((2,nf,nq))
    dp=np.zeros_like(e);power=np.zeros(nf);number=np.zeros(nf)
    for face in range(nf):
        left,right=edges[face,0],edges[face,1]
        low=max(e[left,0],e[right,0]);high=min(e[left,-1],e[right,-1])
        if high<=low:raise ValueError('No shared electronic energy support at a face')
        for k in range(nq):
            energy=(low+high)/2+(high-low)*x[k]/2
            q[face,k]=energy
            measure=(high-low)*gw[k]/2
            for side in range(2):
                node=left if side==0 else right
                a=amplitude[node]
                query=math.sqrt((energy-a)*(energy+a))
                up=np.searchsorted(xi,query,side='right')
                if up<1 or up>=ne:raise ValueError('Transport quadrature lacks native brackets')
                lo=up-1;elo=e[node,lo];ehi=e[node,up]
                b=(energy-elo)/(ehi-elo)
                if b<0 or b>1:raise ValueError('Invalid transport barycentric activity')
                lower[side,face,k]=lo;upper[side,face,k]=up;beta[side,face,k]=b
                lpi[side,face,k]=_interp_log(lp[node],lo,up,b)
                lhi[side,face,k]=_interp_log(lh[node],lo,up,b)
            value=4*face_rates[face]*measure*_difference(
                lpi[0,face,k]+lhi[1,face,k],lpi[1,face,k]+lhi[0,face,k])
            flux[face,k]=value
            for side in range(2):
                node=left if side==0 else right;sign=-1 if side==0 else 1
                lo=lower[side,face,k];up=upper[side,face,k];b=beta[side,face,k]
                dp[node,lo]+=sign*(1-b)*value/(4*w[lo]*volumes[node])
                dp[node,up]+=sign*b*value/(4*w[up]*volumes[node])
            power[face]+=energy*value
            number[face]+=value
        power[face]*=N0*Eu*Eu
        number[face]*=N0*Eu
    return dp,power,number,flux,q,lower,upper,beta,lpi,lhi


@njit(cache=False,inline='always')
def _record(lp,lh,node,k,value):
    if value<0:lp[node,k]-=value
    elif value>0:lh[node,k]+=value


@njit(cache=False,inline='always')
def _fraction(available,loss):
    if loss<=0:return 1.
    q=available/loss
    return min(1.,.95*q)


@njit(cache=False,inline='always')
def _efactor(fp,fh,node,k,value):
    if value<0:return fp[node,k]
    if value>0:return fh[node,k]
    return 1.


@njit(cache=False,parallel=True)
def limit_events(rate,omega,lo,up,beta,i,j,rec,flux,q,tlo,tup,tbeta,
                 edges,vol,w,phcap,node_energy,N0,Eu,oldp,oldn,dt,bp,bh,bn):
    nn,ne=oldp.shape;npop=oldn.shape[1];nev=len(i);nf=len(edges);nq=flux.shape[1]
    lp=np.zeros_like(oldp);lh=np.zeros_like(oldp);ln=np.zeros_like(oldn)
    for node in prange(nn):
        for ev in range(nev):
            value=rate[node,ev]
            if value==0:continue
            _record(lp,lh,node,i[ev],-value if rec[ev] else value)
            _record(lp,lh,node,j[ev],-value)
            if value<0:
                ln[node,lo[node,ev]]-=(1-beta[node,ev])*value
                ln[node,up[node,ev]]-=beta[node,ev]*value
    for face in range(nf):
        for k in range(nq):
            value=flux[face,k]
            for side in range(2):
                node=edges[face,side];sign=-1 if side==0 else 1
                b=tbeta[side,face,k]
                _record(lp,lh,node,tlo[side,face,k],sign*(1-b)*value/vol[node])
                _record(lp,lh,node,tup[side,face,k],sign*b*value/vol[node])
    fp=np.ones_like(oldp);fh=np.ones_like(oldp);fn=np.ones_like(oldn)
    for node in range(nn):
        for k in range(ne):
            cap=4*w[k]
            availablep=bp[node,k]
            availableh=bh[node,k]
            if availablep<0 or availableh<0:
                raise ValueError('Supplied heating exhausts midpoint native-hole budget; reduce dt')
            fp[node,k]=_fraction(availablep,lp[node,k])
            fh[node,k]=_fraction(availableh,lh[node,k])
        for k in range(npop):
            available=bn[node,k]
            if available<0:raise ValueError('Supplied heating exhausts midpoint native-hole budget; reduce dt')
            fn[node,k]=_fraction(available,ln[node,k])
    cdp=np.zeros_like(oldp);cdn=np.zeros_like(oldn);tdp=np.zeros_like(oldp)
    power=np.zeros(nf);number=np.zeros(nf)
    changed=0.;total=0.;count=0;active=0;reduction=0.
    for node in prange(nn):
        for ev in range(nev):
            value=rate[node,ev]
            if value==0:continue
            ii,jj=i[ev],j[ev]
            factor=min(_efactor(fp,fh,node,ii,-value if rec[ev] else value),
                       _efactor(fp,fh,node,jj,-value))
            if value<0:
                if beta[node,ev]<1:factor=min(factor,fn[node,lo[node,ev]])
                if beta[node,ev]>0:factor=min(factor,fn[node,up[node,ev]])
            new=value*factor
            cdp[node,ii]+= -new if rec[ev] else new
            cdp[node,jj]-=new
            cdn[node,lo[node,ev]]+=(1-beta[node,ev])*new
            cdn[node,up[node,ev]]+=beta[node,ev]*new
            unit=node_energy[node]*dt*omega[node,ev]
            changed+=unit*abs(value-new);total+=unit*abs(value)
            if value!=0:
                active+=1
                if factor<1:count+=1
            reduction=max(reduction,1-factor)
        for k in range(ne):cdp[node,k]/=4*w[k]
        for k in range(npop):cdn[node,k]/=phcap[k]
    for face in range(nf):
        for k in range(nq):
            value=flux[face,k];factor=1.
            for side in range(2):
                node=edges[face,side];sign=-1 if side==0 else 1;b=tbeta[side,face,k]
                factor=min(factor,_efactor(fp,fh,node,tlo[side,face,k],sign*(1-b)*value))
                factor=min(factor,_efactor(fp,fh,node,tup[side,face,k],sign*b*value))
            new=value*factor
            for side in range(2):
                node=edges[face,side];sign=-1 if side==0 else 1
                il,iu=tlo[side,face,k],tup[side,face,k];b=tbeta[side,face,k]
                tdp[node,il]+=sign*(1-b)*new/(4*w[il]*vol[node])
                tdp[node,iu]+=sign*b*new/(4*w[iu]*vol[node])
            power[face]+=q[face,k]*new;number[face]+=new
            unit=N0*Eu*Eu*dt*q[face,k]
            changed+=unit*abs(value-new);total+=unit*abs(value)
            if value!=0:
                active+=1
                if factor<1:count+=1
            reduction=max(reduction,1-factor)
        power[face]*=N0*Eu*Eu;number[face]*=N0*Eu
    return cdp,cdn,tdp,power,number,changed,total,count/max(active,1),reduction


"""Node-local storage for the same native acoustic event algebra.
The legacy kernels remain available as a numerical reference.
"""

def ordered_incidence(edges,node_count):
    counts=np.bincount(np.asarray(edges).ravel(),minlength=node_count)
    offsets=np.r_[0,np.cumsum(counts)].astype(np.int64)
    faces=np.empty(len(edges)*2,np.int64);sides=np.empty_like(faces)
    cursor=offsets[:-1].copy()
    for face,(left,right) in enumerate(edges):
        for side,node in enumerate((left,right)):
            index=cursor[node];faces[index]=face;sides[index]=side;cursor[node]+=1
    return offsets,faces,sides

@njit(cache=False,inline='always')
def _record_local(lp,lh,k,value):
    if value<0:lp[k]-=value
    elif value>0:lh[k]+=value

@njit(cache=False,parallel=True)
def acoustic_unlimited(e,amplitude,p,n,lp,lh,ln,l1n,i,j,rec,pairweights,
                       phenergy,phcapacity,w,coefficient,cutlow,cuthigh,prefactor,equilibrium_mask):
    nn,ne=e.shape;npop=n.shape[1];nev=len(i)
    dp=np.zeros_like(p);dn=np.zeros_like(n);invalid=np.zeros(nn,np.int64)
    for node in prange(nn):
        if equilibrium_mask[node]:continue
        a2=amplitude[node]**2
        for ev in range(nev):
            ii,jj=i[ev],j[ev];ei,ej=e[node,ii],e[node,jj]
            om=ei+ej if rec[ev] else ej-ei
            if om<cutlow or om>cuthigh:continue
            up=max(1,np.searchsorted(phenergy,om))
            if up>=npop or om<phenergy[0]:invalid[node]|=1;continue
            lo=up-1;b=(om-phenergy[lo])/(phenergy[up]-phenergy[lo])
            coh=1+a2/(ei*ej) if rec[ev] else 1-a2/(ei*ej)
            c=prefactor*pairweights[ev]*coefficient*om*om*coh
            emit=_interp_log(l1n[node],lo,up,b);absorb=_interp_log(ln[node],lo,up,b)
            if rec[ev]:f=lp[node,ii]+lp[node,jj]+emit;r=lh[node,ii]+lh[node,jj]+absorb
            else:f=lh[node,ii]+lp[node,jj]+emit;r=lp[node,ii]+lh[node,jj]+absorb
            value=c*_difference(f,r)
            dp[node,ii]+= -value if rec[ev] else value
            dp[node,jj]-=value
            dn[node,lo]+=(1-b)*value;dn[node,up]+=b*value
        for k in range(ne):dp[node,k]/=4*w[k]
        for k in range(npop):dn[node,k]/=phcapacity[k]
    if np.any(invalid):raise ValueError('Coupled event lacks phonon bracket; extend physical quadrature')
    return dp,dn

@njit(cache=False,parallel=True)
def acoustic_limited_nodes(e,amplitude,p,n,lp,lh,ln,l1n,i,j,rec,pairweights,
                       phenergy,phcapacity,w,coefficient,cutlow,cuthigh,prefactor,equilibrium_mask,
                       flux,tlo,tup,tbeta,vol,node_energy,dt,bp,bh,bn,offsets,faces,sides):
    nn,ne=e.shape;npop=n.shape[1];nev=len(i);nq=flux.shape[1]
    cdp=np.zeros_like(p);cdn=np.zeros_like(n);fp=np.ones_like(p);fh=np.ones_like(p)
    invalid=np.zeros(nn,np.int64)
    changed=0.;total=0.;count=0;active=0;reduction=0.
    for node in prange(nn):
        rate=np.zeros(nev);omega=np.zeros(nev);beta=np.zeros(nev)
        lower=np.zeros(nev,np.int64);upper=np.ones(nev,np.int64)
        loss_p=np.zeros(ne);loss_h=np.zeros(ne);loss_n=np.zeros(npop)
        a2=amplitude[node]**2
        if not equilibrium_mask[node]:
            for ev in range(nev):
                ii,jj=i[ev],j[ev];ei,ej=e[node,ii],e[node,jj]
                om=ei+ej if rec[ev] else ej-ei;omega[ev]=om
                if om<cutlow or om>cuthigh:continue
                up=max(1,np.searchsorted(phenergy,om))
                if up>=npop or om<phenergy[0]:invalid[node]|=1;continue
                lo=up-1;b=(om-phenergy[lo])/(phenergy[up]-phenergy[lo])
                lower[ev]=lo;upper[ev]=up;beta[ev]=b
                coh=1+a2/(ei*ej) if rec[ev] else 1-a2/(ei*ej)
                c=prefactor*pairweights[ev]*coefficient*om*om*coh
                emit=_interp_log(l1n[node],lo,up,b);absorb=_interp_log(ln[node],lo,up,b)
                if rec[ev]:f=lp[node,ii]+lp[node,jj]+emit;r=lh[node,ii]+lh[node,jj]+absorb
                else:f=lh[node,ii]+lp[node,jj]+emit;r=lp[node,ii]+lh[node,jj]+absorb
                value=c*_difference(f,r);rate[ev]=value
                if value==0:continue
                _record_local(loss_p,loss_h,ii,-value if rec[ev] else value)
                _record_local(loss_p,loss_h,jj,-value)
                if value<0:
                    loss_n[lo]-=(1-b)*value;loss_n[up]-=b*value
        # Exactly the legacy order for each node: collision events, then faces
        # in ascending global face index and quadrature index.
        for incident in range(offsets[node],offsets[node+1]):
            face=faces[incident];side=sides[incident];sign=-1 if side==0 else 1
            for k in range(nq):
                value=flux[face,k];b=tbeta[side,face,k]
                _record_local(loss_p,loss_h,tlo[side,face,k],sign*(1-b)*value/vol[node])
                _record_local(loss_p,loss_h,tup[side,face,k],sign*b*value/vol[node])
        fn=np.ones(npop)
        for k in range(ne):
            availablep=bp[node,k];availableh=bh[node,k]
            if availablep<0 or availableh<0:invalid[node]|=2
            fp[node,k]=_fraction(availablep,loss_p[k]);fh[node,k]=_fraction(availableh,loss_h[k])
        for k in range(npop):
            available=bn[node,k]
            if available<0:invalid[node]|=2
            fn[k]=_fraction(available,loss_n[k])
        for ev in range(nev):
            value=rate[ev]
            if value==0:continue
            ii,jj=i[ev],j[ev]
            factor=min(_efactor(fp,fh,node,ii,-value if rec[ev] else value),
                       _efactor(fp,fh,node,jj,-value))
            if value<0:
                if beta[ev]<1:factor=min(factor,fn[lower[ev]])
                if beta[ev]>0:factor=min(factor,fn[upper[ev]])
            new=value*factor
            cdp[node,ii]+= -new if rec[ev] else new;cdp[node,jj]-=new
            cdn[node,lower[ev]]+=(1-beta[ev])*new;cdn[node,upper[ev]]+=beta[ev]*new
            unit=node_energy[node]*dt*omega[ev]
            changed+=unit*abs(value-new);total+=unit*abs(value)
            if value!=0:
                active+=1
                if factor<1:count+=1
            reduction=max(reduction,1-factor)
        for k in range(ne):cdp[node,k]/=4*w[k]
        for k in range(npop):cdn[node,k]/=phcapacity[k]
    if np.any((invalid&1)!=0):raise ValueError('Coupled event lacks phonon bracket; extend physical quadrature')
    if np.any((invalid&2)!=0):raise ValueError('Supplied heating exhausts midpoint native-hole budget; reduce dt')
    return cdp,cdn,fp,fh,changed,total,count,active,reduction

@njit(cache=False)
def apply_transport(flux,q,tlo,tup,tbeta,edges,vol,w,N0,Eu,dt,fp,fh,changed,total,count,active,reduction):
    nf=len(edges);nq=flux.shape[1]
    tdp=np.zeros_like(fp);power=np.zeros(nf);number=np.zeros(nf)
    for face in range(nf):
        for k in range(nq):
            value=flux[face,k];factor=1.
            for side in range(2):
                node=edges[face,side];sign=-1 if side==0 else 1;b=tbeta[side,face,k]
                factor=min(factor,_efactor(fp,fh,node,tlo[side,face,k],sign*(1-b)*value))
                factor=min(factor,_efactor(fp,fh,node,tup[side,face,k],sign*b*value))
            new=value*factor
            for side in range(2):
                node=edges[face,side];sign=-1 if side==0 else 1
                il,iu=tlo[side,face,k],tup[side,face,k];b=tbeta[side,face,k]
                tdp[node,il]+=sign*(1-b)*new/(4*w[il]*vol[node])
                tdp[node,iu]+=sign*b*new/(4*w[iu]*vol[node])
            power[face]+=q[face,k]*new;number[face]+=new
            unit=N0*Eu*Eu*dt*q[face,k]
            changed+=unit*abs(value-new);total+=unit*abs(value)
            if value!=0:
                active+=1
                if factor<1:count+=1
            reduction=max(reduction,1-factor)
        power[face]*=N0*Eu*Eu;number[face]*=N0*Eu
    return tdp,power,number,changed,total,count/max(active,1),reduction
