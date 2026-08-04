"""
Why does E2 settle? Compute the Jacobian eigenvalues of the (E2, FSH) subsystem
at its fixed point, with the GnRH drive held at its slow-average constant.

If Re(eigenvalues) < 0 for all of them -> the fixed point attracts -> E2 settles.
If they are complex -> it settles by SPIRALLING (damped focus): the "one
overshoot then flat" we saw. The imaginary part gives the damped ring period;
the real part gives the decay time. Neither can be positive here because the
trace is structurally -(kE2+kFSH) < 0.
"""

import numpy as np
from efp_two_loop import fsh_to_e2, e2_to_fsh, kE2, kFSH, cFSH

CE2, K_E2, DRIVE = 0.012, 400.0, 5.0   # from the freeze test

def dfdx(f, x, h=1e-4):
    return (f(x + h) - f(x - h)) / (2 * h)

# --- fixed point of the frozen 2-var system -------------------------------
# dE2  = CE2*fsh_to_e2(FSH) - kE2*E2 = 0  ->  E2 = (CE2/kE2)*fsh_to_e2(FSH)
# dFSH = cFSH*(e2_to_fsh(E2)+DRIVE) - kFSH*FSH = 0
# solve by simple iteration
E2, FSH = 372.0, 7.0
for _ in range(20000):
    E2 = (CE2 / kE2) * fsh_to_e2(FSH)
    FSH = (cFSH / kFSH) * (e2_to_fsh(E2, K_E2) + DRIVE)
print(f"fixed point:  E2* = {E2:.3f}   FSH* = {FSH:.4f}")

# --- Jacobian ------------------------------------------------------------
a11 = -kE2
a12 =  CE2  * dfdx(fsh_to_e2, FSH)
a21 =  cFSH * dfdx(lambda e: e2_to_fsh(e, K_E2), E2)
a22 = -kFSH
J = np.array([[a11, a12], [a21, a22]])

print("\nJacobian at fixed point:")
print(f"  d(dE2)/dE2  = {a11:+.5f}     d(dE2)/dFSH = {a12:+.5f}   (FSH stimulates E2, >0)")
print(f"  d(dFSH)/dE2 = {a21:+.5f}   d(dFSH)/dFSH = {a22:+.5f}   (E2 suppresses FSH, <0)")

tr = np.trace(J); det = np.linalg.det(J)
eig = np.linalg.eigvals(J)
print(f"\ntrace = {tr:+.5f}  (= -(kE2+kFSH), structurally < 0)")
print(f"det   = {det:+.6e}  (> 0)")
print(f"eigenvalues = {eig[0]:.5f} , {eig[1]:.5f}")

re = eig.real; im = eig.imag
print(f"\nRe(eig) max = {re.max():+.5f}  -> {'STABLE (settles)' if re.max()<0 else 'UNSTABLE'}")
if abs(im[0]) > 1e-9:
    period = 2*np.pi/abs(im[0])
    tau_decay = -1/re.max()
    print(f"complex pair -> DAMPED FOCUS (spirals in):")
    print(f"  ring period    ~ {period:8.1f} min ({period/60:.1f} h)")
    print(f"  decay time 1/|Re| ~ {tau_decay:8.1f} min ({tau_decay/60:.1f} h)")
    print(f"  -> overshoots ~{tau_decay/period:.1f} ring(s) worth before flat  = 'one wobble then settle'")
else:
    print("real eigenvalues -> node (settles without overshoot)")
