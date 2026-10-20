"""Topography for the app from the slim MATLAB files (topmodel_explorer_twi_batch.m).

    python scripts/twi_prep.py <folder with twi_slim_<gauge>.mat>     (run from build/)

Writes data/twi.json per gauge:
  ti     30 equal-width TWI classes (as ti_classes() in R/model_utils.R), lnLambdaLin, stats
  route  200-bin histogram of D8 channel distances to the outlet (as routing_distances() in R)
  map    display grid (<= 320 px): per pixel the 0/50/75/90/97/100th percentile TWI class of
         its 30 m cells (uint8, 255 = outside) and a hillshade (uint8), base64
"""
import scipy.io as sio, numpy as np, json, glob, base64, os, sys
from scipy.ndimage import uniform_filter
out={}
CAP=320
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G=json.load(open(os.path.join(ROOT,'scripts','gauges.json'))); IDS=set(G['virginia']+G['other'])
folder=sys.argv[1] if len(sys.argv)>1 else '.'
for f in sorted(glob.glob(os.path.join(folder,'twi_slim_*.mat'))):
    if os.path.basename(f)[9:17] not in IDS: print('skip',f); continue
    M=sio.loadmat(f,squeeze_me=True,struct_as_record=False)
    g=str(M['gauge']); t=M['twiWS'].astype(float); ok=~np.isnan(t); dx=float(M['dx'])
    v=t[ok]
    # --- ti_classes(v, 30) exactly as in model_utils.R
    n=30; e=np.linspace(v.min(), v.max()+1e-9, n+1)
    b=np.clip(np.searchsorted(e, v, side='right'),1,n); b[v>=e[-1]]=n   # findInterval, rightmost closed
    cnt=np.bincount(b,minlength=n+1)[1:]; sm=np.bincount(b,weights=v,minlength=n+1)[1:]
    keep=cnt>0; val=(sm[keep]/cnt[keep]); frac=cnt[keep]/len(v)
    remap=-np.ones(n+1,int); remap[1:][keep]=np.arange(keep.sum())
    mx=v.max(); lnLL=mx+np.log(np.mean(np.exp(v-mx)))
    cls=np.full(t.shape,255,np.uint8); cls[ok]=remap[b]
    # --- routing: hist_counts(distCh[ws], 200)
    d=M['distCh'].astype(float)[ok]; d=d[~np.isnan(d)]
    de=np.linspace(d.min(), d.max()+1e-9, 201); db=np.clip(np.searchsorted(de,d,side='right'),1,200)
    dc=np.bincount(db,minlength=201)[1:]; ctr=(de[1:]+de[:-1])/2
    # --- display grid: blocks of k x k cells
    H,W=t.shape; k=int(np.ceil(max(H,W)/CAP)); h2,w2=int(np.ceil(H/k)),int(np.ceil(W/k))
    P=np.full((h2*k,w2*k),255,np.uint8); P[:H,:W]=cls
    blk=P.reshape(h2,k,w2,k).transpose(0,2,1,3).reshape(h2,w2,k*k).astype(float); blk[blk==255]=np.nan
    inside=np.sum(~np.isnan(blk),axis=2)>=max(1,(k*k)//3)
    qs=[0,50,75,90,97,100]
    Q=np.full((h2,w2,len(qs)),255,np.uint8)
    with np.errstate(all='ignore'):
        qq=np.nanpercentile(blk,qs,axis=2)
    for i in range(len(qs)):
        a=np.where(inside,np.round(qq[i]),255); Q[...,i]=np.nan_to_num(a,nan=255).astype(np.uint8)
    # hillshade from block-mean elevation
    z=M['zWS'].astype(float); Z=np.full((h2*k,w2*k),np.nan); Z[:H,:W]=z
    zb=np.nanmean(Z.reshape(h2,k,w2,k),axis=(1,3))
    zf=np.where(np.isnan(zb),np.nanmean(zb),zb); cs=dx*k
    gy,gx=np.gradient(zf,cs); dzdx=gx; dzdy=-gy
    slope=np.arctan(np.hypot(dzdx,dzdy)); aspect=np.arctan2(-dzdx,-dzdy)
    az,alt=np.radians(315),np.radians(45)
    hs=np.clip(np.sin(alt)*np.cos(slope)+np.cos(alt)*np.sin(slope)*np.cos(az-aspect),0,1)
    hsB=np.where(inside,np.round(hs*254),255).astype(np.uint8)
    out[g]=dict(ti=dict(val=[round(x,4) for x in val],frac=[round(x,6) for x in frac],lnLambdaLin=round(float(lnLL),4),
                    lambda_=round(float(v.mean()),4),sd=round(float(v.std()),3),n=int(len(v)),dx=dx,
                    min=round(float(v.min()),2),max=round(float(v.max()),2)),
                route=dict(dist=[round(x,1) for x in ctr],frac=[round(x,6) for x in dc/dc.sum()],max=round(float(d.max()),0),mean=round(float(d.mean()),0)),
                map=dict(w=w2,h=h2,k=k,px_m=cs,qs=qs,q=base64.b64encode(Q.tobytes()).decode(),hs=base64.b64encode(hsB.tobytes()).decode()),
                area=round(float(len(v)*dx*dx/1e6),2), dem=str(M['stats'].dem_source))
    print(g,'classes',keep.sum(),'lam %.2f'%v.mean(),'map',w2,'x',h2,'k',k,'route max %.0f km'%(d.max()/1e3))
s=json.dumps(out,separators=(',',':'))
open(os.path.join(ROOT,'data','twi.json'),'w').write(s); print('bytes',len(s))
