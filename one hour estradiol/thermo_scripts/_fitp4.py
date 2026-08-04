import numpy as np, pickle
from PIL import Image
from scipy.optimize import minimize
fold=np.load('fold.npy',allow_pickle=True).item()
L=['non-preg','pregnant','Esr1i']
hrs=fold[L[0]][0]; OBS=np.vstack([fold[l][1] for l in L]); AMP=OBS.max(1)-OBS.min(1)
def sig(z): return 1/(1+np.exp(-np.clip(z,-60,60)))
c=np.cos(2*np.pi*(hrs-12.)/24.)
PREG=np.array([0.,1.,1.]); DAYS=[13,14,15,16,17]; mid=np.array(DAYS,float)+0.5
P4M=np.outer(PREG,np.maximum(0.,1.-(mid-13.)/6.5))          # (3,5)
SRC=r'C:\Users\holly\source\garfinkel-xiong-menstrual-cycle\one hour estradiol\mousy temp.png'
img=np.asarray(Image.open(SRC).convert('RGB'))
R_,G_,B_=[img[...,k].astype(float) for k in range(3)]
red=(R_>110)&(R_-G_>45)&(R_-B_>45); red[:250,:]=False
XSPAN=[(116,662),(739,1285),(1362,1907)]
YANCH=[[(326.5,38.0),(365.5,37.0),(403.5,36.0)],[(312.5,38.0),(367.5,37.0),(423.0,36.0)],[(320.5,38.0),(372.0,37.0),(423.0,36.0)]]
omean=[];raw=[]
for i in range(3):
    x0,x1=XSPAN[i]; rows,temps=zip(*YANCH[i]); m,cc=np.polyfit(rows,temps,1)
    xs,ys=[],[]
    for xc in range(x0,x1+1):
        col=np.nonzero(red[:,xc])[0]
        if col.size: xs.append(xc); ys.append(col.mean())
    t=13.0+5.0*(np.array(xs,float)-x0)/(x1-x0); T=m*np.array(ys)+cc
    raw.append((t,T)); omean.append([T[(t>=d)&(t<d+1)].mean() for d in DAYS])
omean=np.array(omean); pickle.dump(raw,open('_raw.pkl','wb')); np.save('_omean.npy',omean)
def waves(z):
    bL,gL,wL,bK,gK,wK,rho,A_T,T0,pL,pK,dr=z
    E2=np.array([0.,1.,rho])[:,None,None]
    p4=P4M[:,:,None]; cb=c[None,None,:]
    dd=(mid-15.5)[None,:,None]
    return T0+A_T*(sig(bK+gK*cb+wK*E2+pK*p4)-sig(bL+gL*cb+wL*E2+pL*p4))+dr*dd
def obj(z):
    W=waves(z); Wm=W.mean(axis=1); r=Wm-OBS
    e=float(np.sqrt(np.mean((np.sqrt(np.mean(r**2,1))/AMP)**2)))
    a=Wm.max(1)-Wm.min(1); ae=float(np.sqrt(np.mean(((a-AMP)/AMP)**2)))
    dm=W.mean(axis=2); de=float(np.sqrt(np.mean(((dm-dm.mean(1,keepdims=True))-(omean-omean.mean(1,keepdims=True)))**2)))
    return e+2.0*ae+1.5*de
B=[(-30,10),(0,60),(0,60),(-30,10),(0,60),(0,60),(0.,1.),(0.3,8),(30,42),(0,25),(0,25),(-0.4,0.1)]
rng=np.random.default_rng(23); lo=np.array([q[0] for q in B]); hi=np.array([q[1] for q in B]); best=None
for _ in range(250):
    z0=lo+rng.random(len(B))*(hi-lo)
    r=minimize(obj,z0,bounds=B,method='L-BFGS-B',options=dict(maxiter=400,ftol=1e-13))
    if best is None or r.fun<best[0]: best=(float(r.fun),r.x)
e,z=best; np.save('best_p4_in_arms.npy',z)
W=waves(z); Wm=W.mean(axis=1); np.save('_W.npy',W)
bL,gL,wL,bK,gK,wK,rho,A_T,T0,pL,pK,dr=z
print(f'pL={pL:.2f} pK={pK:.2f}  drift={dr:+.3f} C/day  rho={rho:.3f}  A_T={A_T:.2f} T0={T0:.2f}')
print(f'gains gL={gL:.1f} gK={gK:.1f} | wL={wL:.1f} wK={wK:.1f}')
print(f'amps model {np.round(Wm.max(1)-Wm.min(1),2)}  observed {np.round(AMP,2)}')
dm=W.mean(axis=2)
for i,l in enumerate(L):
    print(f'  {l:10s} mean drift model {dm[i][-1]-dm[i][0]:+.2f}  observed {omean[i][-1]-omean[i][0]:+.2f}')
