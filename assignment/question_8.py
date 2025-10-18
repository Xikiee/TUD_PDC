import numpy as np
# import pandas as pd

import scipy

import matplotlib
import matplotlib.pyplot as plt

import control as ctrl
# Steady state 1 values

# Power
EP_SS1        = 100.0   # Electrical power (kW)

# Electrolyte Flows (kg/s)
F1_ELEC_SS1   = 2.0     # Electrolyte flow (inlet cathode)
F2_ELEC_SS1   = 2.0     # Electrolyte flow (inlet anode)
F3_ELEC_SS1   = 2.0     # Electrolyte flow (outlet cathode)
F4_ELEC_SS1   = 2.0     # Electrolyte flow (outlet anode)
F7_ELEC_SS1   = 2.0     # Electrolyte flow (outlet V1)
F8_ELEC_SS1   = 2.0     # Electrolyte flow (outlet V2)

# Gas/Species Flows (kg/s)
F3_H2_SS1     = 0.00045 # Hydrogen flow (produced in S1 cathode)
F_OH_SS1      = 0.00773 # Hydroxide flow (consumed/produced in S1)
F4_O2_SS1     = 0.0036  # Oxygen flow (produced in S1 anode)
F5_SS1        = 0.00045  # Hydrogen flow (outlet V1)
F6_SS1        = 0.0036  # Oxygen flow (outlet V2)
F9_SS1        = 0.0     # Water flow (inlet V1)

F5_const = F3_H2_SS1
F6_const = F4_O2_SS1
# Temperatures (°C)
T1_SS1        = 70.0    # Inlet temp. cathode
T2_SS1        = 70.0    # Inlet temp. anode
T3_SS1        = 72.70   # Outlet temp. cathode
T4_SS1        = 72.70   # Outlet temp. anode
TCA_SP_SS1    = 70.0    # Set point TCA

# Time
t_span = [0,int(200*60)]
time = np.linspace(t_span[0],t_span[-1],len(range(t_span[-1]))+1)

EP_array = EP_SS1 * np.ones(len(time))

for i, t in enumerate(time):
    if (20*60) < t <= (30*60):
        EP_array[i] = EP_array[i-1] + 40/60
    elif (90*60) < t <= (100*60):
        EP_array[i] = EP_array[i-1] - 40/60
    else:
        EP_array[i] = EP_array[i-1]



# other constants

Cs1      = 7000.0        # kJ/C
Cp_elec  = 3.24          # kJ/(kg·C)
V1gas    = 0.10          # m^3 (H2)
V2gas    = 0.10          # m^3 (O2)
p2       = 0.002/440.0   # kg/kJ (H2 production per kJ)
p3       = 0.016/440.0   # kg/kJ (O2 production per kJ)
p4       = 0.35          # (fraction of EP to heat)
F1       = F1_ELEC_SS1   # kg/s
F2       = F2_ELEC_SS1   # kg/s

F5_const = F3_H2_SS1
F6_const = F4_O2_SS1

#for ideal gas law
R     = 8.314            # J/(mol·K)
T_gas = 353.0            # K  (≈ 80 C)
M_H2  = 0.002016         # kg/mol
M_O2  = 0.032000         # kg/mol

# 3) Time grid seconds
dt_sec = 1
dt_min = dt_sec/ 60.0         # simulate with the same step in seconds
N = len(time)
t_min = time.copy()/60
t_sec = time.copy()



# Interpolate EP to this grid
EP_series = EP_array.copy()    # kW == kJ/s in our equations


# 4) Allocate & initial conditions
T1    = np.empty(N); T1[0]    = TCA_SP_SS1
T3    = np.empty(N); T3[0]    = T3_SS1
rhoH2 = np.empty(N); rhoH2[0] = 200000/(R*T_gas/M_H2)  #only the first value is of importance, every other value is random
rhoO2 = np.empty(N); rhoO2[0] = 200000/(R*T_gas/M_H2)

print(rhoH2[-1])   # there are move values for densities, and they vary more, idk if this will change stuff later
# For derived plots
F3H2_series = np.empty(N)
PH2 = np.empty(N)
PO2 = np.empty(N)



for k in range(N-1):
    EP_k = EP_series[k]                    # kW (kJ/s)
    # ODEs
    dT1_dt    = (1.0/60.0)*(TCA_SP_SS1 - T1[k])                           # C/s
    Q_stack   = (F1+F2)*Cp_elec*(T1[k] - T3[k]) + p4*EP_k             # kJ/s
    dT3_dt    = Q_stack / Cs1                                          # C/s
    drhoH2_dt = (p2*EP_k - F5_const) / V1gas                           # kg/(m^3*s)
    drhoO2_dt = (p3*EP_k - F6_const) / V2gas                           # kg/(m^3*s)

    # Integrate with dt = 60 s
    T1[k+1]    = T1[k]    + dT1_dt*dt_sec
    T3[k+1]    = T3[k]    + dT3_dt*dt_sec
    rhoH2[k+1] = max(0.0, rhoH2[k] + drhoH2_dt*dt_sec)
    rhoO2[k+1] = max(0.0, rhoO2[k] + drhoO2_dt*dt_sec)

    # Derived at step k (for plotting)
    F3H2_series[k] = p2*EP_k

    # last derived values


F3H2_series[-1] = p2*EP_series[-1]
PH2 = (rhoH2*R*T_gas/M_H2)/1000.0    # kPa
PO2 = (rhoO2*R*T_gas/M_O2)/1000.0    # kPa

print(PH2)       # there are move values for densities, and they vary more, idk if this will change stuff later

# Controller test

I=0
F5_const = F3_H2_SS1

EP_series = EP_array.copy()

gain_H2= -2.89e6
gain_O2= -1.81e5

tau= 30

#rhoH2[0]=

TCA_SP=TCA_SP_SS1

old_level_controller_error=0

for k in range(N-1):
    EP_k = EP_series[k]                    # kW (kJ/s)
    # ODEs
    dT1_dt    = (1.0/60.0)*(TCA_SP - T1[k])                           # C/s
    Q_stack   = (F1+F2)*Cp_elec*(T1[k] - T3[k]) + p4*EP_k             # kJ/s
    dT3_dt    = Q_stack / Cs1                                          # C/s
    drhoH2_dt = (p2*EP_k - F5_const) / V1gas                           # kg/(m^3*s)
    drhoO2_dt = (p3*EP_k - F6_const) / V2gas                           # kg/(m^3*s)

    # Integrate with dt = 60 s
    T1[k+1]    = T1[k]    + dT1_dt*dt_sec
    T3[k+1]    = T3[k]    + dT3_dt*dt_sec
    rhoH2[k+1] = max(0.0, rhoH2[k] + drhoH2_dt*dt_sec)
    rhoO2[k+1] = max(0.0, rhoO2[k] + drhoO2_dt*dt_sec)

    PH2[k+1] = max(0.0, rhoH2[k] + drhoH2_dt*dt_sec)*(R*T_gas/M_H2)/1000.0
    PO2[k+1] = max(0.0, rhoO2[k] + drhoO2_dt*dt_sec)*(R*T_gas/M_O2)/1000.0

    # Derived at step k (for plotting)
    F3H2_series[k] = p2*EP_k

    Kc=2/(gain_H2*tau)
    ti=2*tau

    level_controller_error = (80-T3[k+1])/(1000*3)
    I+=level_controller_error

    if old_level_controller_error * level_controller_error < 0:
      I = 0.0 # Reset the integral term!
      # print(f"I reset at t={k*dt_sec/60}") # Optional: for debugging

    #print(I)
    u_controller = Kc*(level_controller_error + I/ti)


    #print(u_controller)

    TCA_SP+=u_controller
    TCA_SP = max(0.0, TCA_SP)
    old_level_controller_error = level_controller_error
    #print(F5_const)
    # u_out = np.clip(u_controller, 0, 1)
    #u_out = np.clip(u_controller + u0, 0, 1.0)






import matplotlib.pyplot as plt
import numpy as np

# Assuming you have the following data already defined:
# time, EP_series, T1, T3, rhoH2, rhoO2, F5_const, F3H2_series

# Create 2x2 subplots
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.tight_layout(pad=4.0)

# Plot 1: EP
axes[0, 0].plot(time/60, EP_series)
axes[0, 0].set_xlim(0, 200)
axes[0, 0].set_xticks(np.arange(0, 201, 20))
axes[0, 0].set_ylim(0, 600)
axes[0, 0].set_yticks(np.arange(0, 601, 100))
axes[0, 0].set_xlabel("Time (min)")
axes[0, 0].set_ylabel("Power (kW)")
axes[0, 0].set_title("Figure 1: EP profile")
axes[0, 0].grid(True, linestyle='-', alpha=0.6)

# Plot 2: T1 & T3
axes[0, 1].plot(time/60, T1, label="T1")
axes[0, 1].plot(time/60, T3, label="T3")
axes[0, 1].set_xlim(0, 200)
axes[0, 1].set_xticks(np.arange(0, 201, 20))
axes[0, 1].set_ylim(65, 85)
axes[0, 1].set_xlabel("Time (min)")
axes[0, 1].set_ylabel("Temperature (°C)")
axes[0, 1].set_title("Figure 2: T1 & T3 Profiles")
axes[0, 1].legend()
axes[0, 1].grid(True, linestyle='-', alpha=0.6)

# Plot 3: PH2 & PO2
axes[1, 0].plot(time/60, PH2, label="rho_H2")
axes[1, 0].plot(time/60, PO2, label="rho_O2")
axes[1, 0].set_xlim(0, 200)
axes[1, 0].set_xticks(np.arange(0, 201, 20))
axes[1, 0].set_xlabel("Time (min)")
axes[1, 0].set_ylabel("Pressure")
axes[1, 0].set_title("Figure 3: rho_H2 & rho_O2 profiles (from ρ, T=353 K)")
axes[1, 0].legend()
axes[1, 0].grid(True, linestyle='-', alpha=0.6)

# Plot 4: F5 (const) & F3H2
axes[1, 1].plot(time/60, np.full_like(time, F5_const), label="F5 (const)")
axes[1, 1].plot(time/60, F3H2_series, label="F3H2 = p2·EP")
axes[1, 1].set_xlim(0, 200)
axes[1, 1].set_xticks(np.arange(0, 201, 20))
axes[1, 1].set_xlabel("Time (min)")
axes[1, 1].set_ylabel("Flow [kg/s]")
axes[1, 1].set_title("Figure 4: F5 and F3H2 profiles")
axes[1, 1].legend()
axes[1, 1].grid(True, linestyle='-', alpha=0.6)

plt.show()
