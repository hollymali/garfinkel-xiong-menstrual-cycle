"""
Licinio image is in mid-late follicular phase, so GnRH->LH->E2 is used.

    dGnRH/dt = -kG*GnRH   (+ a_n at each pulse)
    dLH/dt   =  cL*GnRH - kL*LH
    dE2/dt   =  cE*LH   - kE*E2
    a_n      = A*[circ] / (1 + (E2d/K_amp)^n_amp)                # NEGATIVE amplitude fb
    period   = Pmin + (Pmax-Pmin)/(1 + (E2d/K_freq)^n_freq)      # POSITIVE frequency fb

"""
# Import numpy
import numpy as np

# First, we define the half-livs
kG = np.log(2) / 4.0     # GnRH t1/2 = 4 min
kL = np.log(2) / 20.0    # LH t1/2 = 20 min
kE = np.log(2) / 30.0    # E2 t/2 = 30 min

# Next, we define the constants
# Solved for so that E2 lands in realistic range
cL = 0.5 
cE = .1145
# Half-saturation densities grounded in E2 mean over mid-late follicular
K = 22.0
A=8.0


# Default run-length is 3 days (for transient to settle) plus 6 days (length of mid-late follicular phase)
def run(n_amp, n_freq, tau=0.0, t_end=12960, dt=0.25, e0=15, period_range=(40,110), circ=False):
    
    steps = int(t_end/dt) # How many steps we will take
    ndelay = int(tau/dt) # How many steps ago time delay is
    time = np.arange(steps+1)*dt # Convert steps back to min

    # Steps + 1 because first spot will be the initial value
    GnRH_array = np.empty(steps + 1) # Create empty array that will store GnRH vals
    LH_array = np.empty(steps + 1) # Same thing for LH
    E2_array = np.empty(steps + 1) # Same thing for E2
    GnRH_array[0] = 0.0 # Initial value of GnRH will be 0 (not on pulse)
    LH_array[0] = 5.0 # Initial value of LH upon entering mid-late follicular phase
    E2_array[0] = e0 # Initial value is specified to test
    phi = 0.0 # Clock that when crossing integers fires a GnRH pulse
    period_min, period_max = period_range

    # Euler's method!
    for i in range(steps):
        delay_index = i - ndelay # Index of tau time ago
        tau_e2 = E2_array[delay_index] if delay_index>=0 else e0 # Grab e2 tau time ago
        period = period_min + (period_max-period_min)/(1 + (tau_e2/K)**n_freq) # NEGATIVE FEEDBACK PERIOD SIGMOID USING TAU E2
        new_phi = phi + dt/period # Approaches integer faster when period is smaller
        fire = np.floor(new_phi) > np.floor(phi) # Evaluates to true when phi crosses whole number
        phi = new_phi # update phi for next step
        if fire:
            amp = A/(1 + (tau_e2/K)**n_amp) # POSITIVE FEEDBACK ON AMPLITUDE
        else:
            amp = 0.0 # If pulse has not fired, no amplitude added
        
        # Get current hormone values
        GnRH = GnRH_array[i]
        LH = LH_array[i]
        E2 = E2_array[i]
        # Now we calculate new value of hormones
        GnRH_array[i+1] = GnRH + dt*(-kG*GnRH) + amp # Add amplitude of pulse and decay to current GnRH
        LH_array[i+1] = LH + dt*(cL*GnRH - kL*LH) # Add inflow and decay to LH
        E2_array[i+1] = E2 + dt*(cE*LH - kE*E2) # Add inflow and decay to E2

    # Return all values (time, hormones)
    return time, GnRH_array, LH_array, E2_array

