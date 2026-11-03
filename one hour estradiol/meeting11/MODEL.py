"""HPG model with a fixed one-hour ovarian response.

Time is in hours from 07:00; SCN phase 16.5 means a 23:30 clock-time peak.
All continuous rates use start-of-step forward Euler updates. GnRH jumps
at integer phase crossings. The optional Windows helper speeds up sliders;
its equations match the Python solver below.
"""

import hashlib
import os
from pathlib import Path
import subprocess
import tempfile
from functools import lru_cache

import numpy as np


REFERENCE = {
    "constant_inflow": 20,
    "production_scale": 250,
    "activation_midpoint": 1,
    "activation_exponent": 4,
    "inhibition_midpoint": 1,
    "inhibition_exponent": 4,
    "n_gain": 1,
    "n_e2_midpoint": 50,
    "central_hill_exponent": 3,
    "n_relaxation": 1,
    "gnrh_pulse_jump": 1,
    "lh_gnrh_gain": 1,
    "period_min": 65,
    "period_max": 75,
    "scn_min": 50,
    "scn_max": 250,
    "sleep_multiplier": 1.57,
    "scn_phase": 16.5,
    "gnrh_half_life_min": 4,
    "lh_half_life_min": 29,
    "e2_half_life_min": 22.5,
}
KEYS = tuple(REFERENCE)
STEP_H = 1.5 / 3600.
WARMUP_H = 96.
OUTPUT_H = 72.


def simulate_endocrine(
    output_h=OUTPUT_H,
    warmup_h=WARMUP_H,
    step_h=STEP_H,
    constant_inflow=20,
    production_scale=250,
    activation_midpoint=1,
    activation_exponent=4,
    inhibition_midpoint=1,
    inhibition_exponent=4,
    n_gain=1,
    n_e2_midpoint=50,
    central_hill_exponent=3,
    n_relaxation=1,
    gnrh_pulse_jump=1,
    lh_gnrh_gain=1,
    period_min=65,
    period_max=75,
    scn_min=50,
    scn_max=250,
    sleep_multiplier=1.57,
    scn_phase=16.5,
    gnrh_half_life_min=4,
    lh_half_life_min=29,
    e2_half_life_min=22.5,
):
    """Integrate N, GnRH, LH, ovarian responsiveness R, E2, and phase."""
    if not (1 <= activation_exponent <= 8 and 1 <= inhibition_exponent <= 8):
        raise ValueError('LH Hill exponents must be in [1,8]')
    if not 0 < step_h <= min(1., 1./n_relaxation,
        gnrh_half_life_min/(60*np.log(2)),lh_half_life_min/(60*np.log(2)),
        e2_half_life_min/(60*np.log(2))):
        raise ValueError('Euler step cannot preserve all state bounds')
    time_h=np.arange(round((warmup_h+output_h)/step_h)+1)*step_h
    n_state = np.empty(len(time_h))
    gnrh = np.empty(len(time_h))
    lh = np.empty(len(time_h))
    e2 = np.empty(len(time_h))
    phi = np.empty(len(time_h))
    responsiveness=np.empty(len(time_h))
    sleep=(((time_h+7)%24>=23)|((time_h+7)%24<7)).astype(float)
    gnrh[0],lh[0],e2[0],phi[0]=0.,3.,15.,0.
    n_state[0]=n_gain*n_e2_midpoint**central_hill_exponent/(n_e2_midpoint**central_hill_exponent+e2[0]**central_hill_exponent)
    responsiveness[0]=inhibition_midpoint**inhibition_exponent/(inhibition_midpoint**inhibition_exponent+lh[0]**inhibition_exponent)
    for i in range(len(time_h)-1):
        phi[i+1]=phi[i]+60.*step_h/(
            (period_max-(period_max-period_min)*n_state[i]**central_hill_exponent/(1.+n_state[i]**central_hill_exponent))
            *(1.+(sleep_multiplier-1.)*sleep[i]))
        n_state[i+1]=n_state[i]+step_h*n_relaxation*(
            n_gain*n_e2_midpoint**central_hill_exponent/(n_e2_midpoint**central_hill_exponent+e2[i]**central_hill_exponent)-n_state[i])
        gnrh[i+1]=gnrh[i]-step_h*(60.*np.log(2.)/gnrh_half_life_min)*gnrh[i]+(
            gnrh_pulse_jump*n_state[i]**central_hill_exponent/(1.+n_state[i]**central_hill_exponent)
            *((scn_max+scn_min)/2.+(scn_max-scn_min)/2.*np.cos(2.*np.pi*(time_h[i]%24.-scn_phase)/24.))
            *(np.floor(phi[i+1])-np.floor(phi[i])))
        lh[i+1]=lh[i]+step_h*(lh_gnrh_gain*gnrh[i]-(60.*np.log(2.)/lh_half_life_min)*lh[i])
        responsiveness[i+1]=responsiveness[i]+step_h*(
            inhibition_midpoint**inhibition_exponent/(inhibition_midpoint**inhibition_exponent+lh[i]**inhibition_exponent)-responsiveness[i])
        e2[i+1]=e2[i]+step_h*(constant_inflow+production_scale
            *lh[i]**activation_exponent/(activation_midpoint**activation_exponent+lh[i]**activation_exponent)
            *responsiveness[i]-(60.*np.log(2.)/e2_half_life_min)*e2[i])
    return dict(time_h=time_h,N=n_state,GnRH=gnrh,LH=lh,E2=e2,phi=phi,
        R=responsiveness,sleep_indicator=sleep)


def simulate_retained(parameters, seconds=1.5, warm=WARMUP_H, output=OUTPUT_H):
    """Return LH/E2 after warmup, using the fast helper when available."""
    executable = cached_runner()
    if executable is None:
        raw = simulate_endocrine(output_h=output, warmup_h=warm,
                                 step_h=seconds / 3600., **parameters)
        first = round(warm * 3600 / seconds)
        result = {key: raw[key][first:] for key in ("LH", "E2")}
    else:
        with tempfile.TemporaryDirectory(prefix="hpg_") as directory:
            output_path = Path(directory) / "trace.bin"
            command = [str(executable), str(output_path), str(warm), str(output), str(seconds),
                       *(str(parameters[key]) for key in KEYS)]
            subprocess.run(command, check=True, capture_output=True,
                           creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            values = np.fromfile(output_path, dtype="<f8").reshape(-1, 2)
        result = dict(LH=values[:, 0], E2=values[:, 1])
    result["time_h"] = (np.arange(len(result["E2"])) + round(warm * 3600 / seconds)) * seconds / 3600
    if not all(np.isfinite(result[key]).all() for key in ("LH", "E2")):
        raise ValueError("Non-finite numerical state for these settings.")
    return result


@lru_cache(maxsize=1)
def cached_runner():
    """Build the optional Windows helper once, outside the project folder."""
    if os.name != "nt":
        return None
    framework = Path(os.environ.get("WINDIR", "C:/Windows")) / "Microsoft.NET"
    compiler = next((framework / family / "v4.0.30319/csc.exe"
                     for family in ("Framework64", "Framework")
                     if (framework / family / "v4.0.30319/csc.exe").is_file()), None)
    if compiler is None:
        return None
    fingerprint = hashlib.sha256(_RUNNER_SOURCE.encode()).hexdigest()[:16]
    cache = Path(tempfile.gettempdir()) / "meeting11_hpg" / fingerprint
    cache.mkdir(parents=True, exist_ok=True)
    executable = cache / "solver.exe"
    if not executable.is_file():
        source = cache / "solver.cs"
        source.write_text(_RUNNER_SOURCE, encoding="utf-8")
        built = subprocess.run([str(compiler), "/nologo", "/optimize+", "/target:exe",
                                f"/out:{executable}", str(source)], capture_output=True,
                               creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        if built.returncode:
            return None
    return executable


# Fast copy of the same Euler equations; one run and LH/E2 output only.
_RUNNER_SOURCE = r"""
using System;
using System.IO;
using System.Globalization;

public class HpgSolver {
    public static void Main(string[] args) {
        double warm=double.Parse(args[1],CultureInfo.InvariantCulture);
        double output=double.Parse(args[2],CultureInfo.InvariantCulture);
        double dt=double.Parse(args[3],CultureInfo.InvariantCulture)/3600.0;
        double[] p=new double[21];
        for(int k=0;k<p.Length;k++) p[k]=double.Parse(args[k+4],CultureInfo.InvariantCulture);
        double constant_inflow=p[0],production_scale=p[1],activation_midpoint=p[2],activation_exponent=p[3],
            inhibition_midpoint=p[4],inhibition_exponent=p[5],n_gain=p[6],n_e2_midpoint=p[7],
            central_hill_exponent=p[8],n_relaxation=p[9],gnrh_pulse_jump=p[10],lh_gnrh_gain=p[11],
            period_min=p[12],period_max=p[13],scn_min=p[14],scn_max=p[15],sleep_multiplier=p[16],
            scn_phase=p[17],gnrh_half_life_min=p[18],lh_half_life_min=p[19],e2_half_life_min=p[20];
        if(dt>1 || activation_exponent<1 || activation_exponent>8 || inhibition_exponent<1 || inhibition_exponent>8)
            throw new Exception("Invalid state step or LH powers");
        int total=(int)Math.Round((warm+output)/dt),start=(int)Math.Round(warm/dt);
        double N=n_gain*Math.Pow(n_e2_midpoint,central_hill_exponent)/(Math.Pow(n_e2_midpoint,central_hill_exponent)+Math.Pow(15,central_hill_exponent));
        double G=0,L=3,E=15,phi=0,R=Math.Pow(inhibition_midpoint,inhibition_exponent)/(Math.Pow(inhibition_midpoint,inhibition_exponent)+Math.Pow(3,inhibition_exponent));
        double ke_pow=Math.Pow(n_e2_midpoint,central_hill_exponent),
            ka_pow=Math.Pow(activation_midpoint,activation_exponent),
            ki_pow=Math.Pow(inhibition_midpoint,inhibition_exponent),
            kg=60*Math.Log(2)/gnrh_half_life_min,
            kl=60*Math.Log(2)/lh_half_life_min,
            kel=60*Math.Log(2)/e2_half_life_min;
        using(var writer=new BinaryWriter(new BufferedStream(File.Create(args[0]),1048576))) {
            for(int i=0;i<=total;i++){
                if(i>=start){writer.Write(L);writer.Write(E);}
                if(i==total)break;
                double t=i*dt,clock=(t+7)%24,sleep=(clock>=23 || clock<7)?1:0;
                double nh_pow=Math.Pow(N,central_hill_exponent),
                    e_pow=Math.Pow(E,central_hill_exponent),
                    la_pow=Math.Pow(L,activation_exponent),
                    li_pow=Math.Pow(L,inhibition_exponent);
                double phi_next=phi+60*dt/((period_max-(period_max-period_min)*nh_pow/(1+nh_pow))*(1+(sleep_multiplier-1)*sleep));
                double N_next=N+dt*n_relaxation*(n_gain*ke_pow/(ke_pow+e_pow)-N);
                double G_next=G-dt*(kg)*G+gnrh_pulse_jump*nh_pow/(1+nh_pow)*((scn_max+scn_min)/2+(scn_max-scn_min)/2*Math.Cos(2*Math.PI*(t%24-scn_phase)/24))*(Math.Floor(phi_next)-Math.Floor(phi));
                double L_next=L+dt*(lh_gnrh_gain*G-(kl)*L);
                double R_next=R+dt*(ki_pow/(ki_pow+li_pow)-R);
                double E_next=E+dt*(constant_inflow+production_scale*la_pow/(ka_pow+la_pow)*R-(kel)*E);
                N=N_next;G=G_next;L=L_next;E=E_next;R=R_next;phi=phi_next;
                if(double.IsNaN(E)||double.IsInfinity(E)||E<0||L<0||N<0||R<0||R>1)throw new Exception("State bounds failed");
            }
        }
    }
}
"""


if __name__ == "__main__":
    trace = simulate_retained(REFERENCE, output=24.)
    for key, units in (("LH", "IU/L"), ("E2", "pg/mL")):
        values = trace[key]
        print(f"{key} first/mean/last: "
              f"{values[0]:.3f}/{values.mean():.3f}/{values[-1]:.3f} {units}")
