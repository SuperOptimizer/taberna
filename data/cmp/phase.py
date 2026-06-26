#!/usr/bin/env python3
"""phase.py VOL.f32 [SCALE=4]  — CENTER-FREE phase-coherence scorer vs held-out VC3D refs.

Why a new metric: METRIC A (angle(grad W, n)) is LOCAL — a field can have a small local
angle yet a consistently-signed error that ACCUMULATES into large winding drift over a long
sheet. That accumulated drift is exactly the "long coherent surface" failure, and A can't see
it. This scorer measures it directly, with NO umbilicus/center (the previous phase attempts were
confounded by a guessed center).

Key geometry (verified): each VC3D segment is a (u,v) grid of world (z,y,x). The v-direction
(rows) is the scroll AXIS/height (dP/dv ~ pure +z); moving along v stays on the SAME physical
wrap. So a sheet-tracking field's W must stay ~CONSTANT along v over the whole segment height.
Both u and v are tangent directions (perp to the true normal n), so locally W shouldn't change
in EITHER — but along v there is also no winding accumulation, making it the clean phase axis.

METRIC C — PHASE DRIFT vs step: for grid steps k along v, take pairs (u,v),(u,v+k) both valid
and inside the volume, report |dW| in WRAP UNITS (W is wraps). A perfect field -> 0 at every k.
Drift growing with k = the iso-surface peeling off the sheet over distance. The number is
literally "after tracing k grid-rows (~k*dz_world voxels) along one sheet, the level set has
drifted this many wraps off it." Center-free, monotone-in-badness, directly the stated goal.

Also reports the same along u (winding direction) as a reference: u-drift includes legitimate
winding accumulation where the sheet curves, so it's NOT a pure phase test, but a field that
tracks should still keep LOCAL (small-k) u-drift low.

archive_lod coord = (world-origin)/SCALE-[z0,y0,x0].  origin = ROI trim (256,1280,768) Paris4.
"""
import numpy as np, tifffile as tf, glob, sys, os

VOL=sys.argv[1]; SCALE=int(sys.argv[2]) if len(sys.argv)>2 and sys.argv[2].isdigit() else 4
OZ,OY,OX=256,1280,768

with open(VOL,'rb') as f:
    dz,dy,dx,z0,y0,x0=[int(v) for v in np.fromfile(f,np.int32,6)]
    W=np.fromfile(f,np.float32,dz*dy*dx).reshape(dz,dy,dx)
fin=np.isfinite(W)

segs=sorted(glob.glob('/home/forrest/paris4_segments/auto_grown_*'))

def sampleW(P):
    """P[...,3] world zyx -> (W at rounded local idx, inside-mask). Material-only (finite)."""
    loc=(P-np.array([OZ,OY,OX]))/SCALE-np.array([z0,y0,x0])
    li=np.round(loc).astype(int)
    ins=(li[...,0]>=0)&(li[...,0]<dz)&(li[...,1]>=0)&(li[...,1]<dy)&(li[...,2]>=0)&(li[...,2]<dx)
    out=np.full(P.shape[:-1],np.nan,np.float64)
    z=np.clip(li[...,0],0,dz-1); y=np.clip(li[...,1],0,dy-1); x=np.clip(li[...,2],0,dx-1)
    w=W[z,y,x]; ok=ins&np.isfinite(w)
    out[ok]=w[ok]
    return out

STEPS=[1,2,4,8,16,32,64]
v_drift={k:[] for k in STEPS}; u_drift={k:[] for k in STEPS}
for d in segs:
    X=tf.imread(d+'/x.tif');Y=tf.imread(d+'/y.tif');Z=tf.imread(d+'/z.tif')
    Vm=(X>=0)&np.isfinite(X)
    P=np.stack([Z,Y,X],-1).astype(np.float64)
    Ws=sampleW(P)                      # W sampled on the whole grid (nan where outside/air)
    valid=Vm&np.isfinite(Ws)
    H,Wd=X.shape
    for k in STEPS:
        # along v (rows, axis 0) = height/same-wrap : the phase test
        if H>k:
            a=Ws[:H-k]; b=Ws[k:]; m=valid[:H-k]&valid[k:]
            dvv=np.abs(b-a)[m]
            if dvv.size: v_drift[k].append(dvv)
        # along u (cols, axis 1) = winding direction : reference
        if Wd>k:
            a=Ws[:,:Wd-k]; b=Ws[:,k:]; m=valid[:,:Wd-k]&valid[:,k:]
            duu=np.abs(b-a)[m]
            if duu.size: u_drift[k].append(duu)

print(f"\n== PHASE {os.path.basename(VOL)}  (SCALE={SCALE}, region z{z0}+{dz} y{y0}+{dy} x{x0}+{dx}) ==")
print("  step |  v-drift (SAME-WRAP axis; want ~0, flat in k) |  u-drift (winding axis; ref)")
print("       |   median    p90     n                          |  median    p90")
for k in STEPS:
    v=np.concatenate(v_drift[k]) if v_drift[k] else np.array([])
    u=np.concatenate(u_drift[k]) if u_drift[k] else np.array([])
    vs=f"{np.median(v):.3f}   {np.percentile(v,90):.3f}  {v.size:>8}" if v.size else "   -"
    us=f"{np.median(u):.3f}   {np.percentile(u,90):.3f}" if u.size else "   -"
    print(f"  {k:>4} |  {vs}   |  {us}")
# headline: phase drift at the largest reliable step, in wrap units
kbig=max(k for k in STEPS if v_drift[k])
vb=np.concatenate(v_drift[kbig])
print(f"\n  HEADLINE phase-drift @ v-step {kbig}: median={np.median(vb):.3f} wraps, "
      f"p90={np.percentile(vb,90):.3f} wraps  (0 = level set never leaves the sheet)")
