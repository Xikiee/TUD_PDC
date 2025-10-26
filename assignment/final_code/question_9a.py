import numpy as np
import pandas as pd

import scipy

import matplotlib
import matplotlib.pyplot as plt
import control as ctrl
# STEP TEST

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
F3_H2_SS1     = 100*0.002/440.0 # Hydrogen flow (produced in S1 cathode)
F_OH_SS1      = 0.00773 # Hydroxide flow (consumed/produced in S1)
F4_O2_SS1     = 100*0.016/440.0  # Oxygen flow (produced in S1 anode)
F5_SS1        = 100*0.002/440.0  # Hydrogen flow (outlet V1)
F6_SS1        = 100*0.016/440.0 # Oxygen flow (outlet V2)
F9_SS1        = 0.0     # Water flow (inlet V1)

# Temperatures (°C)
T1_SS1        = 70.0    # Inlet temp. cathode
T2_SS1        = 70.0    # Inlet temp. anode
T3_SS1        = 72.70   # Outlet temp. cathode
T4_SS1        = 72.70   # Outlet temp. anode
TCA_SP_SS1    = 70.0    # Set point TCA

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
dt_min = 1.0
dt_sec = dt_min * 60.0         # simulate with the same step in seconds
N = len(times)                 # 201 points
t_min = times.copy()
t_sec = t_min * 60.0

# Interpolate EP to this grid
EP_series = EP_array.copy()    # kW == kJ/s in our equations


# 4) Allocate & initial conditions
T1    = np.empty(N); T1[0]    = TCA_SP_SS1
T3    = np.empty(N); T3[0]    = T3_SS1
rhoH2 = np.empty(N); rhoH2[0] = 200000/(R*T_gas/M_H2)
rhoO2 = np.empty(N); rhoO2[0] = 200000/(R*T_gas/M_H2)

# For derived plots
F3H2_series = np.empty(N)
PH2 = np.empty(N)
PO2 = np.empty(N)


for k in range(N-1):
    EP_k = 100.0                   # kW (kJ/s)

    if k==19:
      print(t_min[k])
      #F5_const = F3_H2_SS1*0.95
      #F6_const = F4_O2_SS1*0.95
      TCA_SP_SS1 = 70.0*1.05

    # ODEs
    dT1_dt    = (1.0/60.0)*(TCA_SP_SS1 - T1[k])                           # C/s
    Q_stack   = (F1+F2)*Cp_elec*(T1[k] - T3[k]) + p4*EP_k             # kJ/s
    dT3_dt    = Q_stack / Cs1                                          # C/s
    drhoH2_dt = (p2*EP_k - F5_const) / V1gas                           # kg/(m^3*s)
    drhoO2_dt = (p3*EP_k - F6_const) / V2gas                           # kg/(m^3*s)

    #print(p3*EP_k-F6_const)

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

#print(PH2)



# Plot 2: T1 & T3
plt.figure(figsize=(7,4))
plt.plot(t_min, T1, label="T1")
plt.plot(t_min, T3, label="T3")
plt.axvline(19, color = 'k', linestyle = '--', label = 'Step at 150 min')
plt.xlim(0, 200); plt.xticks(np.arange(0, 201, 20))
plt.xlabel("Time (min)"); plt.ylabel("Temperature (C0)"); plt.title(" T3 Profile")
plt.legend(); plt.grid(True, linestyle='-', alpha=0.6)

print(T3[20])
print(f'gain:',((T3[199]-T3[10])/(70*0.05)))


# Calculate the 63.2% target value
T3_initial = T3[10] # Using T3[10] as your initial steady state
T3_final = T3[199]
y_63_2_target = T3_initial + (T3_final - T3_initial) * 0.632


i_tau = np.argmin(np.abs(T3 - y_63_2_target))

# The time constant (tau) estimation would then be:
t_63_2 = t_min[i_tau]

# --- Example Output Check ---
print(f"Target T3 Value (Y_63.2%): {y_63_2_target}")
print(f"Index of Closest Match (i_tau): {i_tau}")
print(f"T3 Value at that Index: {T3[i_tau]}")
print(f"Time (t_63.2) at that Index: {t_63_2}")

# Plot 3: PH2 & PO2
plt.figure(figsize=(7,4))
plt.plot(t_min, rhoH2, label="rho_H2")
plt.plot(t_min, rhoO2, label="rho_O2")
plt.xlim(0, 200); plt.xticks(np.arange(0, 201, 20))
plt.xlabel("Time [min]"); plt.ylabel("denisty"); plt.title("Figure 3: rho_H2 & rho_O2 profiles (from ρ, T=353 K)")
plt.legend(); plt.grid(True, linestyle='-', alpha=0.6)