#!/usr/bin/env python3
"""wrap.py VOL.f32 [SCALE=4]  — CENTER-FREE wrap-counting / phase scorer vs held-out VC3D refs.

The decisive test for "long coherent surfaces": does the field assign a CONSISTENT winding
increment between two physically-adjacent loops of a real sheet? No umbilicus needed — the GT
geometry itself defines "one wrap":

  Each VC3D segment is a single physical sheet that spirals ~10x around the scroll. For a point p
  on it (with true normal n), the radially-ADJACENT loop of the SAME segment sits ~pitch away
  ALONG n (tangentially close, normal-offset ~ one sheet spacing). Those two points are exactly
  ONE WRAP apart by construction. A field that counts wraps cleanly must give W(q)-W(p) = a
  consistent value (ideally 1) at EVERY such pair. Spread / non-integer / sign-flips = the field
  miscounts = leaks & switches = bad long surfaces.

Outputs (all center-free):
  pitch_world  : measured sheet spacing (median normal-offset to adjacent loop) — set dr_per_winding≈pitch/SCALE
  |dW| median  : winding increment across one real wrap (want ~1.0; <1 undercount/leak, >1 overcount)
  frac in[.7,1.3]: fraction counting one wrap correctly
  CV           : spread of dW (lower = more uniform counting)
  integer-conc : circular concentration R of frac(dW) (1=every pair an integer apart=coherent, 0=random)
  sign-consistency: fraction of pairs whose dW sign agrees with the majority (1=no fold mis-orderings)

archive_lod coord = (world-origin)/SCALE-[z0,y0,x0].  origin = ROI trim (256,1280,768) Paris4.
"""
import numpy as np, tifffile as tf, glob, sys, os
from scipy.spatial import cKDTree

VOL=sys.argv[1]; SCALE=int(sys.argv[2]) if len(sys.argv)>2 and sys.argv[2].isdigit() else 4
OZ,OY,OX=256,1280,768
with open(VOL,'rb') as f:
    dz,dy,dx,z0,y0,x0=[int(v) for v in np.fromfile(f,np.int32,6)]
    W=np.fromfile(f,np.float32,dz*dy*dx).reshape(dz,dy,dx)

def sampleW(P):
    loc=(P-np.array([OZ,OY,OX]))/SCALE-np.array([z0,y0,x0]); li=np.round(loc).astype(int)
    ins=(li[...,0]>=0)&(li[...,0]<dz)&(li[...,1]>=0)&(li[...,1]<dy)&(li[...,2]>=0)&(li[...,2]<dx)
    z=np.clip(li[...,0],0,dz-1);y=np.clip(li[...,1],0,dy-1);x=np.clip(li[...,2],0,dx-1)
    w=W[z,y,x]; out=np.where(ins&np.isfinite(w),w,np.nan); return out

segs=sorted(glob.glob('/home/forrest/paris4_segments/auto_grown_*'))
rng=np.random.default_rng(0)
dWs=[]; pitches=[]
for d in segs:
    X=tf.imread(d+'/x.tif');Y=tf.imread(d+'/y.tif');Z=tf.imread(d+'/z.tif')
    Vm=(X>=0)&np.isfinite(X); P=np.stack([Z,Y,X],-1).astype(np.float64)
    H,Wd=X.shape
    du=P[1:-1,2:]-P[1:-1,:-2]; dv=P[2:,1:-1]-P[:-2,1:-1]
    ok=Vm[1:-1,1:-1]&Vm[1:-1,2:]&Vm[1:-1,:-2]&Vm[2:,1:-1]&Vm[:-2,1:-1]
    n=np.cross(du,dv); ln=np.linalg.norm(n,axis=-1); ok&=ln>1e-6
    n=n/np.where(ln[...,None]>0,ln[...,None],1); C=P[1:-1,1:-1]
    ys,xs=np.where(ok)
    if ys.size<50: continue
    pts=C[ys,xs]; nrm=n[ys,xs]; Wp=sampleW(pts)
    tree=cKDTree(pts)
    nsel=min(4000,pts.shape[0]); sel=rng.choice(pts.shape[0],nsel,replace=False)
    for i in sel:
        if not np.isfinite(Wp[i]): continue
        p=pts[i]; nn=nrm[i]
        idx=tree.query_ball_point(p,60.0)
        if len(idx)<5: continue
        idx=np.array(idx); q=pts[idx]-p
        nc=q@nn; tc=np.linalg.norm(q-np.outer(nc,nn),axis=1)
        m=(tc<15)&(nc>8)                # adjacent loop on the +n side only (consistent handedness
        if not m.any(): continue        # within a segment) -> dW sign is a real fold/mis-order test
        jj=np.argmin(nc[m]); j=idx[m][jj]   # nearest such = immediately-adjacent +n loop
        if not np.isfinite(Wp[j]): continue
        dw=Wp[j]-Wp[i]
        dWs.append(dw); pitches.append((pts[j]-p)@nn)
dW=np.array(dWs); pit=np.array(pitches)
print(f"\n== WRAP {os.path.basename(VOL)}  (SCALE={SCALE}, region z{z0}+{dz} y{y0}+{dy} x{x0}+{dx})  n={dW.size} adjacent-loop pairs ==")
if dW.size<50: print("  too few pairs (volume z-range misses the GT)"); sys.exit()
a=np.abs(dW); med=np.median(a)
frac=np.mean((a>0.7)&(a<1.3)); cv=a.std()/(a.mean()+1e-9)
fr=dW-np.round(dW); R=np.sqrt(np.mean(np.cos(2*np.pi*fr))**2+np.mean(np.sin(2*np.pi*fr))**2)
signcon=max(np.mean(dW>0),np.mean(dW<0))   # oriented by +n: a coherent field rises consistently
print(f"  pitch_world  median={np.median(pit):.1f}  (=> ideal dr_per_winding ~ {np.median(pit)/SCALE:.1f} vox at SCALE={SCALE})")
print(f"  |dW| per real wrap   median={med:.3f}   (want ~1.0; <1 undercount/leak, >1 overcount)")
print(f"  counts-one-wrap  frac dW in [0.7,1.3] = {100*frac:.0f}%")
print(f"  uniformity       CV(|dW|) = {cv:.2f}    (lower=better)")
print(f"  integer-coherence R(frac dW) = {R:.3f}   (1=every pair an integer apart, 0=random phase)")
print(f"  sign-consistency = {100*signcon:.0f}%    (1=no fold mis-orderings)")
