"""How much of the temperature model is actually earning its keep?

Fits three models to the same folded CBT data with the same loss, and compares
against the CEILING for any cosine-shaped model (a free cosine per arm).

  CEILING  9 par : CBT_i(t) = M_i - A_i*cos(2pi(t-ph_i)/24)     -- pure description
  SIMPLE   8 par : A_i = A0/(1 + E2_i/K),  M_i = T0 + m*P4_i    -- E2 compresses, P4 offsets
  DUAL    12 par : two opposed sigmoids + first-order plant     -- the current model
"""
import io
import numpy as np
from scipy.optimize import minimize

fold = np.load("fold.npy", allow_pickle=True).item()
LABELS = ["non-preg", "pregnant", "Esr1i"]
hrs = fold[LABELS[0]][0]
OBS = np.vstack([fold[l][1] for l in LABELS])
AMP = OBS.max(axis=1) - OBS.min(axis=1)
NB = OBS.shape[1]

# pull the DUAL model's definitions without running freephase's own fit
src = io.open("freephase.py", encoding="utf-8").read().split("ALL = np.array")[0]
ns = {}
exec(src, ns)
model_dual, MODELS, TAIL = ns["model"], ns["MODELS"], ns["TAIL"]


def loss_from(M):
    r = M - OBS
    return float(np.sqrt(np.mean((np.sqrt(np.mean(r ** 2, axis=1)) / AMP) ** 2)))


# ---------- CEILING: free cosine per arm ----------
def ceiling_model(x):
    M, A, ph = x[0:3], x[3:6], x[6:9]
    return np.vstack([M[i] - A[i] * np.cos(2 * np.pi * (hrs - ph[i]) / 24.0)
                      for i in range(3)])


# ---------- SIMPLE ----------
def simple_model(x):
    T0, m, A0, K, rho = x[:5]
    ph = x[5:8]
    E2 = np.array([0.0, 1.0, rho])
    P4 = np.array([0.0, 1.0, 1.0])
    A = A0 / (1.0 + E2 / max(K, 1e-6))
    M = T0 + m * P4
    return np.vstack([M[i] - A[i] * np.cos(2 * np.pi * (hrs - ph[i]) / 24.0)
                      for i in range(3)])


def fit(fn, bnds, nstart, seed=3):
    rng = np.random.default_rng(seed)
    lo = np.array([b[0] for b in bnds]); hi = np.array([b[1] for b in bnds])
    best = None
    for _ in range(nstart):
        x0 = lo + rng.random(len(bnds)) * (hi - lo)
        r = minimize(lambda z: loss_from(fn(z)), x0, bounds=bnds,
                     method="L-BFGS-B", options=dict(maxiter=400, ftol=1e-14))
        if best is None or r.fun < best[0]:
            best = (float(r.fun), r.x)
    return best


B_CEIL = [(35, 39)] * 3 + [(0, 3)] * 3 + [(0, 24)] * 3
B_SIMP = [(35, 39), (-2, 2), (0, 4), (1e-3, 20), (0, 1)] + [(0, 24)] * 3

e_ceil, x_ceil = fit(ceiling_model, B_CEIL, 300)
e_simp, x_simp = fit(simple_model, B_SIMP, 300)

spec = MODELS["dual       "]
bnds = spec["bounds"] + TAIL
rng = np.random.default_rng(31)
lo = np.array([b[0] for b in bnds]); hi = np.array([b[1] for b in bnds])
ALL = np.array([0, 1, 2])
best = None
for _ in range(300):
    x0 = lo + rng.random(len(bnds)) * (hi - lo)
    r = minimize(lambda z: loss_from(model_dual(spec, z, ALL)), x0, bounds=bnds,
                 method="L-BFGS-B", options=dict(maxiter=400, ftol=1e-14))
    if best is None or r.fun < best[0]:
        best = (float(r.fun), r.x)
e_dual, x_dual = best

print(f"{'model':10s} {'params':>7} {'loss':>8}   amplitudes (obs 2.10 0.72 1.16)")
for nm, e, M, npar in (("CEILING", e_ceil, ceiling_model(x_ceil), 9),
                       ("SIMPLE", e_simp, simple_model(x_simp), 8),
                       ("DUAL", e_dual, model_dual(spec, x_dual, ALL), 12)):
    a = M.max(axis=1) - M.min(axis=1)
    print(f"{nm:10s} {npar:7d} {e:8.4f}   {np.round(a, 2)}")

T0, m, A0, K, rho = x_simp[:5]
print(f"\nSIMPLE fitted: T0={T0:.2f} C   m(P4 offset)={m:+.2f} C   A0={A0:.2f} C   "
      f"K={K:.3f}   rho={rho:.3f}")
print(f"  predicted amplitude ratios: "
      f"{np.round(1.0/(1.0+np.array([0,1,rho])/K) / (1.0/(1.0+0/K)), 2)}   (observed 1.00 0.34 0.55)")
print(f"  phases (h): {np.round(x_simp[5:8], 1)}   (observed peaks 8.2 15.8 13.2)")
