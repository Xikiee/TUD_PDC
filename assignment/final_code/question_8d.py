import numpy as np
import matplotlib.pyplot as plt

# ---- CONSTANTS ----
EP_SS1 = 100  # electrical power 

# Electrolyte Flows (kg/s)
F1_ELEC_SS1 = 2.0     # Electrolyte flow (inlet cathode)
F2_ELEC_SS1 = 2.0     # Electrolyte flow (inlet anode)

# Gas/Species Flows (kg/s)
F3_H2_SS1 = 0.00045  # Hydrogen flow (produced in S1 cathode)
F4_O2_SS1 = 0.0036   # Oxygen flow (produced in S1 anode)
F5_const = F3_H2_SS1
F6_const = F4_O2_SS1

# Temperatures (°C)
T1_SS1 = 70.0   # Inlet temp. cathode
T3_SS1 = 72.70  # Outlet temp. cathode
TCA_SP_SS1 = 70.0  # Set point TCA

# ---- CONSTANT PARAMETERS ----
Cs1 = 7000.0        # kJ/C
Cp_elec = 3.24      # kJ/(kg·C)
p2 = 0.002 / 440.0  # kg/kJ (H2 production per kJ)
p3 = 0.016 / 440.0  # kg/kJ (O2 production per kJ)
p4 = 0.35           # (fraction of EP to heat)

# ---- defining time frame ----
dt = 0.01  # s
t_span = [0, 200 * 60]
time = np.arange(t_span[0], t_span[-1] + dt, dt)
iter = len(time)

# ---- defining EP Array ----
EP_array = EP_SS1 * np.ones(len(time))
for i, t in enumerate(time):
    if (20 * 60) < t <= (30 * 60):
        EP_array[i] = EP_array[i - 1] + 40 / 60 * dt
    elif (90 * 60) < t <= (100 * 60):
        EP_array[i] = EP_array[i - 1] - 40 / 60 * dt
    else:
        EP_array[i] = EP_array[i - 1]
EP = EP_array.copy()

# ---- Define the gains ----
gain_H2 = -3.3e7
gain_O2 = -2.08e6
gain_vec = np.array([gain_H2, gain_O2])

# ---- Define the Disturbance ----
F3 = p2 * EP
F4 = p3 * EP
disturbance = np.array([F3, F4])

# ---- Initial conditions ----
int_flow = np.array([F5_const, F6_const])
pressure0 = [2e5, 2e5]
SP = np.array([2e5, 3e5])

# ---- Set point ---- 
# Pressure: 
def pressure_set(time:float,SP:float)->np.ndarray:
    """
    This function defines the step change of the pressure in time
    Arguments: 
    Time (float): the time frame in which the pressure is to be evaluated
    SP (float): Set point between which the set point will occure. 
set flow is
    Returns: The return is the pressure step change in np.ndarray format 
    """
    pressure = []
    for i in time:
        if i < 150*60:
            pressure.append(SP[0]) 
        else: 
            pressure.append(SP[1])
    return pressure


tau = 30

# ---- Temperature Parameters ----
temp_params = [TCA_SP_SS1, T1_SS1, T3_SS1, F1_ELEC_SS1, F2_ELEC_SS1]

# ---- INTEGRATION FUNCTION ----
def func_int_pressure(
    pressure0:np.ndarray, gain:np.ndarray, disturbance:np.ndarray,
    int_flow:np.ndarray, SP:np.ndarray, dt:float, tau:float,
    temp_params:np.ndarray, iter:int)->np.ndarray:
    
    TCA_SP_SS1, T1_SS1, T3_SS1, F1_ELEC_SS1, F2_ELEC_SS1 = temp_params

    n = len(gain)
    pressure = np.zeros((n, iter))
    pressure[:, 0] = pressure0
    flow = np.zeros_like(pressure)
    flow[:, 0] = int_flow


    # Initialize temperature arrays
    T1 = np.zeros(iter)
    T3 = np.zeros(iter)
    T1[0] = T1_SS1
    T3[0] = T3_SS1

    int_error = np.zeros(n)

    for i in range(1, iter):
        if time[i] < 150 * 60:
            set_pressure = SP[0]
        else:
            set_pressure = SP[1]

        # ---- Pressure dynamics ----
        dP_dt = gain * (flow[:, i - 1] - disturbance[:, i - 1])     # F3 - F5 
        pressure[:, i] = pressure[:, i - 1] + dP_dt * dt

        # ---- Temperature dynamics ----
        dT1_dt = (1 / 60) * (TCA_SP_SS1 - T1[i - 1])
        dT3_dt = (1 / Cs1) * (
            (F1_ELEC_SS1 + F2_ELEC_SS1) * Cp_elec * (T1[i - 1] - T3[i - 1])
            + p4 * EP[i - 1]
        )
        T1[i] = T1[i - 1] + dT1_dt * dt
        T3[i] = T3[i - 1] + dT3_dt * dt

        # ---- Controller ----
        error = set_pressure - pressure[:, i]
        int_error += error
        K_c = 2 / (gain * tau)
        T_c = 2 * tau
        flow[:, i] = flow[:, 0] + error * K_c + (K_c * int_error * dt) / T_c
        # use flow[:,0] because the equation is 

    return pressure, flow, K_c, T_c, T1, T3, int_error


# ---- Run simulation ----
pressure, flow, K_c, T_c, T1, T3, error = func_int_pressure(
    pressure0, gain_vec, disturbance, int_flow, SP, dt, tau, temp_params, iter
)

# ---- Plotting ----
fig, axes = plt.subplots(3, 2, figsize=(10, 6))

# EP profile
axes[0, 0].plot(time / 60, EP_array, label="EP")
axes[0, 0].set_xlabel("Time (min)")
axes[0, 0].set_ylabel("Power (kW)")
axes[0, 0].set_title("EP profile")
axes[0, 0].legend()
axes[0, 0].set_xlim(0, 200)
axes[0, 0].grid(which="both")

# ---- Temperature (T1,T3)
axes[0, 1].plot(time / 60, T1, label="T1")
axes[0, 1].plot(time / 60, T3, label="T3")
axes[0, 1].set_xlabel("Time (min)")
axes[0, 1].set_ylabel("Temperature (°C)")
axes[0, 1].set_title("Temperature vs Time")
axes[0, 1].legend()
axes[0, 1].set_xlim(0, 200)
axes[0, 1].grid(which="both")

# ---- Pressure PC1.SP, PC1.Pv
axes[1, 0].plot(time / 60, pressure[0], label="PH2")
axes[1, 0].plot(time / 60, pressure_set(time, SP), label="Step Change", ls='--', color='orange')
axes[1, 0].set_xlabel("Time (min)")
axes[1, 0].set_ylabel("Pressure (Pa)")
axes[1, 0].set_title("Pressure vs Time: PC1.SP, PC1.PV")
axes[1, 0].legend()
axes[1, 0].set_xlim(0, 200)
axes[1, 0].grid(which="both")

# ---- Pressure PC2.SP, PC2.PV
axes[1, 1].plot(time / 60, pressure[1], label="PO2")
axes[1, 1].plot(time / 60, pressure_set(time, SP), label="Step Change", ls="--", color='orange')
axes[1, 1].set_xlabel("Time (min)")
axes[1, 1].set_ylabel("Pressure (Pa)")
axes[1, 1].set_title("Pressure vs Time: PC2.SP, PC2.PV")
axes[1, 1].legend()
axes[1, 1].set_xlim(0, 200)
axes[1, 1].grid(which="both")

# ---- Flow F5, F3_H2
axes[2, 0].plot(time / 60, flow[0], label="Flow H2 (F5)")
axes[2, 0].plot(time / 60, F3, label="Flow F3", color='orange')
axes[2, 0].set_xlabel("Time (min)")
axes[2, 0].set_ylabel("Flow (kg/s)")
axes[2, 0].set_title("Flow vs Time: F5 vs F3_H2")
axes[2, 0].legend()
axes[2, 0].set_xlim(0, 200)
axes[2, 0].grid(which="both")

# ---- Flow F6, F4_O2
axes[2, 1].plot(time / 60, flow[1], label="Flow O2 (F6)")
axes[2, 1].plot(time / 60, F4, label="Flow F4", color='orange')
axes[2, 1].set_xlabel("Time (min)")
axes[2, 1].set_ylabel("Flow (kg/s)")
axes[2, 1].set_title("Flow vs Time: F6 vs F4_O2")
axes[2, 1].legend()
axes[2, 1].set_xlim(0, 200)
axes[2, 1].grid(which="both")

# --- Compact layout ---
plt.tight_layout(h_pad=0.8)
plt.subplots_adjust(hspace=0.3)
plt.show()

print(f"K_c : {K_c}, and T_c : {T_c}")
