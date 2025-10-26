import numpy as np
import matplotlib.pyplot as plt

# STEP TEST (density) - SECOND-BASED


EP_SS1        = 100.0   # Electrical power (kW)

# Electrolyte Flows (kg/s)
F1_ELEC_SS1   = 2.0
F2_ELEC_SS1   = 2.0
F3_ELEC_SS1   = 2.0
F4_ELEC_SS1   = 2.0
F7_ELEC_SS1   = 2.0
F8_ELEC_SS1   = 2.0

# Gas/Species Flows (kg/s)
F3_H2_SS1     = 100*0.002/440.0
F_OH_SS1      = 0.00773
F4_O2_SS1     = 100*0.016/440.0
F5_SS1        = 100*0.002/440.0
F6_SS1        = 100*0.016/440.0
F9_SS1        = 0.0

# Temperatures (°C)
T1_SS1        = 70.0
T2_SS1        = 70.0
T3_SS1        = 72.70
T4_SS1        = 72.70
TCA_SP_SS1    = 70.0

# Stack & gas properties
Cs1      = 7000.0
Cp_elec  = 3.24
V1gas    = 0.10
V2gas    = 0.10
p2       = 0.002/440.0
p3       = 0.016/440.0
p4       = 0.35
F1       = F1_ELEC_SS1
F2       = F2_ELEC_SS1
F5_const = F3_H2_SS1
F6_const = F4_O2_SS1

R     = 8.314
T_gas = 353.0
M_H2  = 0.002016
M_O2  = 0.032000



t_sec = np.arange(0, 200*60 + 1, 1)  # 0..12000 s
N = len(t_sec)
dt_sec = 1.0
t_min = t_sec / 60.0


EP_series = np.full(N, EP_SS1)  # constant 100 kW for step test


rhoH2 = np.zeros(N); rhoH2[0] = 200000/(R*T_gas/M_H2)
rhoO2 = np.zeros(N); rhoO2[0] = 200000/(R*T_gas/M_O2)
F3H2_series = np.zeros(N)
PH2 = np.zeros(N)
PO2 = np.zeros(N)


for k in range(N-1):
    EP_k = EP_series[k]  # kW

    # Step perturbation at t = 150 min
    if t_min[k] == 150:
        F5_const *= 1.01
        F6_const *= 1.01

    # ODEs
    drhoH2_dt = (p2*EP_k - F5_const) / V1gas
    drhoO2_dt = (p3*EP_k - F6_const) / V2gas

    # Integrate (1-second step)
    rhoH2[k+1] = max(0.0, rhoH2[k] + drhoH2_dt * dt_sec)
    rhoO2[k+1] = max(0.0, rhoO2[k] + drhoO2_dt * dt_sec)

    # Derived
    F3H2_series[k] = p2*EP_k

# Last derived value
F3H2_series[-1] = p2*EP_series[-1]

# Calculate pressures
PH2 = rhoH2 * R * T_gas / M_H2 / 1000.0  # kPa
PO2 = rhoO2 * R * T_gas / M_O2 / 1000.0  # kPa


idx_150 = 150*60
idx_200 = 200*60
step = 0.01

slope_h2 = (rhoH2[idx_200] - rhoH2[idx_150]) / (t_sec[idx_200] - t_sec[idx_150])
slope_O2 = (rhoO2[idx_200] - rhoO2[idx_150]) / (t_sec[idx_200] - t_sec[idx_150])

Gain_h2 = slope_h2 / (step * F5_SS1)
Gain_O2 = slope_O2 / (step * F6_SS1)

print("Gain H2:", Gain_h2)
print("Gain O2:", Gain_O2)
print("Slope H2:", slope_h2)
print("Slope O2:", slope_O2)

fig, ax = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

# Left plot: H2
ax[0].plot(t_min, rhoH2, label="rho_H2", color='red')
ax[0].axvline(x=150, color='black', linestyle='--', linewidth=1.5, label='Step at t = 150 min')
ax[0].text(165, 0.10, 'slope = -0.00272', color='black', fontsize=10)
ax[0].set_xlim(0, 200)
ax[0].set_xticks(np.arange(0, 201, 20))
ax[0].set_xlabel("Time (min)")
ax[0].set_ylabel("Density (kg/m^3)")
ax[0].set_title("Density_H2 Profile (1-s simulation)")
ax[0].legend()
ax[0].grid(True, linestyle='-', alpha=0.6)

# Right plot: O2
ax[1].plot(t_min, rhoO2, label="rho_O2", color='blue')
ax[1].axvline(x=150, color='black', linestyle='--', linewidth=1.5, label='Step at t = 150 min')
ax[1].text(165, 2, 'slope = -0.0218', color='black', fontsize=10)
ax[1].set_xlim(0, 200)
ax[1].set_xticks(np.arange(0, 201, 20))
ax[1].set_xlabel("Time (min)")
ax[1].set_title("Density_O2 Profile (1-s simulation)")
ax[1].legend()
ax[1].grid(True, linestyle='-', alpha=0.6)

plt.tight_layout()
plt.show()
