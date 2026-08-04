import numpy as np
from scipy.optimize import minimize
from scipy.signal import find_peaks
fold=np.load('fold.npy',allow_pickle=True).item()
L=['non-preg','pregnant','Esr1i']
hrs=fold[L[0]][0]; OBS=np.vstack([fold[l][1] for l in L]); AMP=OBS.max(1)-OBS.min(1)
def sig(z): return 1/(1+np.exp(-np.clip(z,-60,60)))
c=np.cos(2*np.pi*(hrs-12.)/24.)
E2=np.array([0.,1.,0.])[:,None]      # Esr1i transduces NO estrogen
P4=np.array([0.,1.,1.])[:,None]      # both pregnant animals have progesterone
NB=len(hrs)
def model(z):
    bL,gL,wL,bK,gK,wK,pL,pK,A_T,T0=z
    return T0+A_T*(sig(bK+gK*c+wK*E2+pK*P4)-sig(bL+gL*c+wL*E2+pL*P4))
def obj(z):
    M=model(z); r=M-OBS
    e=float(np.mean(np.sqrt(np.mean(r**2,1))/AMP))
    sd=float(np.mean(((M.std(1)-OBS.std(1))/OBS.std(1))**2))
    edge=sum(1.0 for i in range(3) for k in (int(np.argmax(M[i])),int(np.argmin(M[i]))) if k<2 or k>NB-3)
    return e+2.0*sd+1.0*edge
B=[(-30,10),(0,60),(0,60),(-30,10),(0,60),(0,60),(0,25),(0,25),(0.3,8),(30,42)]
rng=np.random.default_rng(77); lo=np.array([q[0] for q in B]); hi=np.array([q[1] for q in B]); best=None
for _ in range(400):
    z0=lo+rng.random(len(B))*(hi-lo)
    r=minimize(obj,z0,bounds=B,method='L-BFGS-B',options=dict(maxiter=500,ftol=1e-13))
    if best is None or r.fun<best[0]: best=(float(r.fun),r.x)
e,z=best; np.save('best_noE2_esr.npy',z); M=model(z)
bL,gL,wL,bK,gK,wK,pL,pK,A_T,T0=z
print(f'bL={bL:.2f} gL={gL:.1f} wL={wL:.1f} pL={pL:.2f} | bK={bK:.2f} gK={gK:.1f} wK={wK:.1f} pK={pK:.2f} | A_T={A_T:.2f} T0={T0:.2f}')
print()
print(f'{"animal":10s}{"peaks":>26}{"amp mod":>9}{"amp obs":>9}{"SD mod":>8}{"SD obs":>8}')
for i,l in enumerate(L):
    a=np.ptp(M[i]); ext=np.concatenate([M[i]]*3)
    pk,_=find_peaks(ext,prominence=0.15*max(a,1e-9))
    pt=[hrs[p-NB] for p in pk if NB<=p<2*NB]
    print(f'{l:10s}{str(np.round(pt,1)):>26}{a:9.2f}{AMP[i]:9.2f}{M[i].std():8.3f}{OBS[i].std():8.3f}')
print()
print('OBSERVED peaks: non-preg 8.2,15.2 | pregnant 15.8,21.8 | Esr1i 13.2 ONLY')
