"""Quick visual diagnostics for menstrual_model.py (verification only)."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import menstrual_model as m


def panel_fig(phase, sol, gnrh, title, fname, zoom=600):
    t = sol.t / 60.0  # hours
    series = [('GnRH', gnrh), ('FSH', sol.y[1]), ('LH', sol.y[2]),
              ('E2', sol.y[0]), ('P4', sol.y[3])]
    fig, ax = plt.subplots(5, 2, figsize=(13, 11), sharex='col')
    fig.suptitle(title, fontsize=13)
    nz = sol.t < zoom
    for r, (nm, y) in enumerate(series):
        ax[r, 0].plot(t, y, lw=0.6)
        ax[r, 0].set_ylabel(nm)
        ax[r, 1].plot(t[nz], y[nz], lw=0.8)
    ax[0, 0].set_title('full run (5.5 days)')
    ax[0, 1].set_title(f'zoom: first {zoom} min')
    ax[-1, 0].set_xlabel('time (hours)')
    ax[-1, 1].set_xlabel('time (hours)')
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(fname, dpi=110)
    plt.close(fig)
    print('wrote', fname)


# 1) EFP faithful (as printed, cFSH = 2e-7)
sol = m.run('EFP')
panel_fig('EFP', sol, m.gnrh_trace('EFP', sol),
          'EFP -- faithful transcription (cFSH = 2e-7 as printed)', 'fig_efp_faithful.png')

# 2) EFP recalibrated (cFSH bumped so FSH enters the feedback-active range)
p = dict(m.EFP); p['cFSH'] = 5e-3
rhs = lambda t, y: m.efp_rhs(t, y, p)
solr = solve_ivp(rhs, (0, 8000), [200, 5, 10, 5, 0], method='RK45',
                 t_eval=np.arange(0, 8000, 0.5), max_step=0.25, rtol=1e-6, atol=1e-9)
amp = np.array([m.efp_amplitude(e) for e in solr.y[0]])
gr = amp * m.pulse_kernel(solr.y[4])
panel_fig('EFP', solr, gr, 'EFP -- recalibrated (cFSH = 5e-3)', 'fig_efp_recal.png')

# 3) MLP faithful
solm = m.run('MLP')
panel_fig('MLP', solm, m.gnrh_trace('MLP', solm),
          'MLP -- faithful transcription', 'fig_mlp.png')
