# Simulation with PI pressure controllers for PC1 (H2) and PC2 (O2)
# This code extends the user's model to add two PI controllers that manipulate F5 and F6
# to keep PH2 and PO2 at 2 bar (200 kPa). A small set-point change is introduced at t=150 min.
# The code produces the requested plots and computes simple performance metrics (IAE, max overshoot).

import numpy as np
import matplotlib.pyplot as plt

# --- User / model parameters (kept from user's script) ---
EP_SS1        = 100.0   # kW
F3_H2_SS1     = 0.00045 # kg/s  (H2 production)
F4_O2_SS1     = 0.0036  # kg/s  (O2 production)
F5_const0     = F3_H2_SS1
F6_const0     = F4_O2_SS1

Cs1      = 7000.0        # kJ/C
Cp_elec  = 3.24          # kJ/(kg·C)
V1gas    = 0.10          # m^3 (H2)
V2gas    = 0.10          # m^3 (O2)
p2       = 0.002/440.0   # kg/kJ (H2 production per kJ)
p3       = 0.016/440.0   # kg/kJ (O2 production per kJ)
p4       = 0.35          # fraction of EP to heat

TCA_SP_SS1    = 70.0
T3_SS1        = 72.70

# ideal gas
R     = 8.314            # J/(mol·K)
T_gas = 353.0            # K  (≈ 80 C)
M_H2  = 0.002016         # kg/mol
M_O2  = 0.032000         # kg/mol

# Time grid (seconds)
dt_sec = 1
sim_min = 200
N = sim_min*60 + 1
time = np.linspace(0, sim_min*60, N)
time_min = time/60.0


EP_array = EP_SS1 * np.ones(len(time))

for i, t in enumerate(time):
    if (20*60) < t <= (30*60):
        EP_array[i] = EP_array[i-1] + 40/60
    elif (90*60) < t <= (100*60):
        EP_array[i] = EP_array[i-1] - 40/60
    else:
        EP_array[i] = EP_array[i-1]

EP = EP_array
# Initial conditions
T1 = np.zeros(N); T1[0] = TCA_SP_SS1
T3 = np.zeros(N); T3[0] = T3_SS1

rhoH2 = np.zeros(N)
rhoO2 = np.zeros(N)
# pick reasonable initial densities to match roughly 2 bar (200 kPa)
# p = rho * R * T / M / 1000  (kPa)
target_p_kPa = 200.0  # 2 bar
rhoH2[0] = target_p_kPa*1000.0 * M_H2 / (R*T_gas)
rhoO2[0] = target_p_kPa*1000.0 * M_O2 / (R*T_gas)

PH2 = np.zeros(N); PO2 = np.zeros(N)
PH2[0] = (rhoH2[0]*R*T_gas/M_H2)/1000.0
PO2[0] = (rhoO2[0]*R*T_gas/M_O2)/1000.0

F3H2_series = np.zeros(N)
F4O2_series = np.zeros(N)

# --- Controllers: choose process gains and tuning parameters ---
# Use the user's approximate gains (pressure per flow) from their script
# Note: these are negative because increasing outlet flow reduces pressure.
gain_H2 = -2.89e6   # [kPa per (kg/s)] approximate (user provided)
gain_O2 = -1.81e5   # [kPa per (kg/s)] approximate (user provided)

# IMC-like tuning used in user's derivation: Kc = 2/(K*tau), Ti = 2*tau  (K is process gain)
tau_c = 30.0  # chosen tuning parameter (seconds) - user used 30 (but units were min? we keep seconds here)
Kc_H2 = 2.0/(gain_H2 * tau_c)
Kc_O2 = 2.0/(gain_O2 * tau_c)
Ti_H2 = 2.0 * tau_c
Ti_O2 = 2.0 * tau_c

# For implementation with 1 s timestep, convert integrator factor accordingly
# We'll implement discrete PI: u = Kc*(e + (1/Ti)*integral_e)
int_H2 = 0.0
int_O2 = 0.0

# Manipulated variables initial (start from steady outlet flows)
F5 = np.ones(N)*F5_const0
F6 = np.ones(N)*F6_const0

# Setpoints (kPa). Both controllers set to 2 bar initially.
SP_H2 = np.ones(N)*200.0  # kPa
SP_O2 = np.ones(N)*200.0  # kPa
# Introduce small SP changes at t = 150 min (9000 s)
sp_change_time = 150*60
for i,t in enumerate(time):
    if t >= sp_change_time:
        SP_H2[i] = 202.0  # +2 kPa step
        SP_O2[i] = 198.0  # -2 kPa step

# Limits on manipulated variables (physical)
F5_min, F5_max = 0.0, 0.01
F6_min, F6_max = 0.0, 0.02

# Storage for controller outputs for plotting
u5 = np.zeros(N)
u6 = np.zeros(N)

# For performance metrics
abs_error_H2 = np.zeros(N)
abs_error_O2 = np.zeros(N)
overshoot_H2 = 0.0
overshoot_O2 = 0.0

# Simulation loop
for k in range(N-1):
    EPk = EP[k]  # kW ~ kJ/s
    # thermal ODEs (unchanged)
    dT1_dt = (1.0/60.0)*(TCA_SP_SS1 - T1[k])  # C/s
    Q_stack = (2.0+2.0)*Cp_elec*(T1[k] - T3[k]) + p4*EPk
    dT3_dt = Q_stack / Cs1
    T1[k+1] = T1[k] + dT1_dt*dt_sec
    T3[k+1] = T3[k] + dT3_dt*dt_sec

    # production flows from EP at time k
    F3H2 = p2 * EPk
    F4O2 = p3 * EPk
    F3H2_series[k] = F3H2
    F4O2_series[k] = F4O2

    # simple mass balance for gases (rho dynamics)
    drhoH2_dt = (F3H2 - F5[k]) / V1gas
    drhoO2_dt = (F4O2 - F6[k]) / V2gas

    rhoH2[k+1] = max(0.0, rhoH2[k] + drhoH2_dt * dt_sec)
    rhoO2[k+1] = max(0.0, rhoO2[k] + drhoO2_dt * dt_sec)

    # measured pressures (kPa)
    PH2[k+1] = (rhoH2[k+1]*R*T_gas/M_H2)/1000.0
    PO2[k+1] = (rhoO2[k+1]*R*T_gas/M_O2)/1000.0

    # Controller calculations each second
    # Error = setpoint - measurement (kPa)
    e_H2 = SP_H2[k+1] - PH2[k+1]
    e_O2 = SP_O2[k+1] - PO2[k+1]
    abs_error_H2[k+1] = abs(e_H2)
    abs_error_O2[k+1] = abs(e_O2)

    # Anti-windup: if sign changes we reset integral (simple approach)
    if k>0 and ( (SP_H2[k]-PH2[k]) * e_H2 < 0 ):
        int_H2 = 0.0
    if k>0 and ( (SP_O2[k]-PO2[k]) * e_O2 < 0 ):
        int_O2 = 0.0

    int_H2 += e_H2 * dt_sec
    int_O2 += e_O2 * dt_sec

    # PI outputs in terms of delta-flow (kg/s). Note sign: increasing outlet flow reduces pressure,
    # so to raise pressure (e>0) controller should reduce outlet flow (negative delta).
    u_H2 = Kc_H2 * (e_H2 + (1.0/Ti_H2) * int_H2)
    u_O2 = Kc_O2 * (e_O2 + (1.0/Ti_O2) * int_O2)

    # compute new manipulated flows
    F5[k+1] = np.clip(F5[k] + u_H2, F5_min, F5_max)
    F6[k+1] = np.clip(F6[k] + u_O2, F6_min, F6_max)

    u5[k+1] = u_H2
    u6[k+1] = u_O2

# Finalize series tails
F3H2_series[-1] = p2*EP[-1]
F4O2_series[-1] = p3*EP[-1]

# Performance metrics
IAE_H2 = np.trapz(abs_error_H2, time)  # Integral Absolute Error (kPa*s)
IAE_O2 = np.trapz(abs_error_O2, time)

max_overshoot_H2 = (np.max(PH2) - 200.0) if np.max(PH2)>200.0 else 0.0
max_overshoot_O2 = (np.max(PO2) - 200.0) if np.max(PO2)>200.0 else 0.0

# --- Plotting similar to report's Table 3 requirement ---
fig, axs = plt.subplots(3,2, figsize=(14,12))
plt.subplots_adjust(hspace=0.4)

# EP profile
axs[0,0].plot(time_min, EP)
axs[0,0].set_title("EP profile")
axs[0,0].set_xlim(0, sim_min)
axs[0,0].set_ylabel("Power (kW)")
axs[0,0].grid(True)

# Temperatures T1, T3
axs[0,1].plot(time_min, T1, label='T1'); axs[0,1].plot(time_min, T3, label='T3')
axs[0,1].set_title("Temperatures")
axs[0,1].legend(); axs[0,1].set_xlim(0, sim_min); axs[0,1].set_ylabel("°C"); axs[0,1].grid(True)

# Pressures and setpoints
axs[1,0].plot(time_min, PH2, label='PH2 (kPa)'); axs[1,0].plot(time_min, SP_H2, '--', label='SP_H2')
axs[1,0].plot(time_min, PO2, label='PO2 (kPa)'); axs[1,0].plot(time_min, SP_O2, '--', label='SP_O2')
axs[1,0].set_title("Pressures and Setpoints")
axs[1,0].legend(); axs[1,0].set_xlim(0, sim_min); axs[1,0].set_ylabel("kPa"); axs[1,0].grid(True)

# Flows F5 and F6
axs[1,1].plot(time_min, F5, label='F5 (kg/s)'); axs[1,1].plot(time_min, F6, label='F6 (kg/s)')
axs[1,1].set_title("Manipulated outlet flows (F5, F6)")
axs[1,1].legend(); axs[1,1].set_xlim(0, sim_min); axs[1,1].set_ylabel("kg/s"); axs[1,1].grid(True)

# Controller outputs (u deltas)
axs[2,0].plot(time_min, u5, label='u5 delta F5'); axs[2,0].plot(time_min, u6, label='u6 delta F6')
axs[2,0].set_title("Controller outputs (delta flows)")
axs[2,0].legend(); axs[2,0].set_xlim(0, sim_min); axs[2,0].set_ylabel("kg/s"); axs[2,0].grid(True)

# Production flows
axs[2,1].plot(time_min, F3H2_series, label='F3H2 (prod)'); axs[2,1].plot(time_min, F4O2_series, label='F4O2 (prod)')
axs[2,1].set_title("Production flows (from EP)")
axs[2,1].legend(); axs[2,1].set_xlim(0, sim_min); axs[2,1].set_ylabel("kg/s"); axs[2,1].grid(True)

plt.show()

# Print performance summary
print("Controller tuning and performance summary")
print(f"Kc_H2 = {Kc_H2:.3e}, Ti_H2 = {Ti_H2:.1f} s")
print(f"Kc_O2 = {Kc_O2:.3e}, Ti_O2 = {Ti_O2:.1f} s")
print(f"IAE_H2 = {IAE_H2:.2f} kPa·s, IAE_O2 = {IAE_O2:.2f} kPa·s")
print(f"Max overshoot PH2 above 200 kPa = {max_overshoot_H2:.3f} kPa")
print(f"Max overshoot PO2 above 200 kPa = {max_overshoot_O2:.3f} kPa")

# Save arrays to disk for further analysis if needed (not required, but convenient)
import os
outdir = '/mnt/data/sim_outputs'
os.makedirs(outdir, exist_ok=True)
np.save(outdir + '/time_min.npy', time_min)
np.save(outdir + '/PH2.npy', PH2)
np.save(outdir + '/PO2.npy', PO2)
np.save(outdir + '/F5.npy', F5)
np.save(outdir + '/F6.npy', F6)

print(f"Saved numpy arrays to {outdir}/")


 