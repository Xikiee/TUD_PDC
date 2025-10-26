import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import scipy 
import control as ctrl

# ---- CONSTANTS ----
EP_SS1 = 100 # electrical power 

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

Cs1      = 7000.0        # kJ/C
Cp_elec  = 3.24          # kJ/(kg·C)
V1gas    = 0.10          # m^3 (H2)
V2gas    = 0.10          # m^3 (O2)
p2       = 0.002/440.0   # kg/kJ (H2 production per kJ)
p3       = 0.016/440.0   # kg/kJ (O2 production per kJ)
p4       = 0.35          # (fraction of EP to heat)
F1       = F1_ELEC_SS1   # kg/s
F2       = F2_ELEC_SS1   # kg/s



# ---- defining time frame ----
dt = 0.01               # smaller dt to have a more accurate integration 
t_span = [0, 200*60]
time = np.arange(t_span[0], t_span[-1] + dt, dt)
iter = len(time)

# ---- defining EP Array ----

EP_array =  EP_SS1* np.ones(len(time))

for i, t in enumerate(time):
    if (20*60) < t <= (30*60):
        EP_array[i] = EP_array[i-1] + 40/60 *dt
    elif (90*60) < t <= (100*60):
        EP_array[i] = EP_array[i-1] - 40/60*dt
    else:
        EP_array[i] = EP_array[i-1]

EP = EP_array.copy()



# ---- Define the gains ----
gain_H2 = -3.3e7
gain_O2 = -2.08e6
gain_vec = np.array([gain_H2, gain_O2])

# ---- Define the Disturbance ----
F3 = p2 * EP
F5 = p3 * EP
disturbance = np.array([F3,F5])

# ---- Define the manipulated variable : the flow ----
F5 = 0.00045  # Hydrogen flow (outlet V1)
F6 = 0.0036  # Oxygen flow (outlet V2)
# manipulated_variable = np.array([F5,F6])

int_flow = np.array([F5,F6])

# ---- Initial pressure ----
pressure0 = [2e5,2e5]

# ---- Set Pressure ----
SP = np.array([2e5,3e5])

# ---- Controller ----
tau = 20




# ---- INTEGRATION FUNCTION : TO INTEGRATE FOR PRESSURE ----
def func_int_pressure(pressure0: np.ndarray, gain: np.ndarray,
                      disturbance: np.ndarray,
                      int_flow:np.ndarray,SP:np.ndarray, dt: float, 
                      tau:float, iter: int) -> np.ndarray:
    

    n = len(gain)
    pressure = np.zeros((n, iter))
    pressure[:, 0] = pressure0

    flow = np.zeros_like(pressure)
    flow[:,0] = int_flow

    int_error = np.zeros(n)
    for i in range(1, iter):

        if time[i] < 150*60:
            set_pressure = SP[0]
        else:
            set_pressure = SP[1]
            
        dP_dt = gain * (flow[:,i-1] - disturbance[:, i-1])  
        pressure[:, i] = pressure[:, i-1] + dP_dt * dt

        # ---- DEFINING THE CONTROLLER ----
        error = set_pressure - pressure[:,i]             
        int_error += error

        K_c = 2/(gain*tau)
        T_c = 2*tau
        flow[:,i]  = flow[:,0] + error*K_c + (K_c*int_error*dt)/T_c                 

    return pressure, flow, K_c, T_c



pressure, flow, K_c, T_c = func_int_pressure(pressure0, gain_vec, disturbance, int_flow,SP,dt,tau,iter)
# print(f"the pressures are {pressure}, the flows are {flow}")


# ---- FIXED PLOTTING CODE ----
fig, axes = plt.subplots(2, 2, figsize=(10, 4))  # create 1x2 grid of subplots

# ---- EP profile ----
axes[0,0].plot(time/60, EP_array, label='EP')
axes[0,0].set_xlabel('Time (min)')
axes[0,0].set_ylabel('Power (kW)')
axes[0,0].set_title('EP profile')
axes[0,0].legend()
axes[0,0].minorticks_on()
axes[0,0].grid(which = 'both')


# --- Pressure plot ---
axes[1,0].plot(time/60, pressure[0], label='PH2')
axes[1,0].plot(time/60, pressure[1], label='PO2')
axes[1,0].set_xlabel('Time (min)')
axes[1,0].set_ylabel('Pressure')
axes[1,0].set_title('Pressure vs Time')
axes[1,0].legend()
axes[1,0].minorticks_on()
axes[1,0].grid(which = 'both')

# --- Flow plot ---
axes[0,1].plot(time/60, flow[0], label='Flow H2')
axes[0,1].plot(time/60, flow[1], label='Flow O2')
axes[0,1].set_xlabel('Time (min)')
axes[0,1].set_ylabel('Flow')
axes[0,1].set_title('Flow vs Time')
axes[0,1].legend()
axes[0,1].minorticks_on()
axes[0,1].grid(which= 'both')

# ---- Temperature plot ----
axes[1,1].plot(time/60, flow[0], label='Flow H2')
axes[1,1].plot(time/60, flow[1], label='Flow O2')
axes[1,1].set_xlabel('Time (min)')
axes[1,1].set_ylabel('Flow')
axes[1,1].set_title('Flow vs Time')
axes[1,1].legend()
axes[1,1].minorticks_on()
axes[1,1].grid(which= 'both')





plt.tight_layout()
plt.show()


print(f"K_c : {K_c}, and T_c : {T_c}")




