import numpy as np
import matplotlib.pyplot as plt
import control as ctrl 

import scipy

####################################### parameters: 
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

# Define the EP array
EP_array=EP_SS1*np.ones(len(time))

for i,t in enumerate(time):
  if (20*60)<t<=(30*60):
    EP_array[i]=EP_array[i-1]+40*60
  elif (90*60)<t<=(100*60):
    EP_array[i]=EP_array[i-1]-(40*60)
  else:
    EP_array[i]=EP_array[i-1]

# Plotting the EP profile
plt.figure(figsize=(10, 5))
plt.plot(time/60,EP_array)
plt.xticks(np.arange(0, t_span[-1]/60, 10))

plt.ylim(0, 600)
plt.yticks(np.arange(0, 601, 100))

plt.xlabel('Time (min)')
plt.ylabel('Power (kW)')
plt.title('Figure 1: The EP profile.')
plt.grid(True, linestyle='-', alpha=0.6)

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
dt_min = 1.0
dt_sec = dt_min * 60.0         # simulate with the same step in seconds
N = len(time)                 # 201 points
t_min = time.copy()
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


########### Question 8 c: IMC controller tuning

#plotting set input at t = 150 min
u = np.zeros_like(time)
u[time >= 150*60] = 1.0

# using the normal f(s) formula 1/(lambda*s+1)**n
# lam = 40 # (lambda 30 looks good)
# num = [1]
# den = [lam**2,2*lam,3]

# for system with pure integration:   
lam = 40
num = [(2*lam), 0, 1]
den = [(lam**2 + 4*lam), 2*lam, 3]
func = ctrl.TransferFunction(num, den)

# Compute response
t, y = ctrl.forced_response(func, T=time, U=u)

# Plot in minutes
plt.figure(figsize=(10, 5))
plt.plot(t/60, y, label='System Response')
plt.xlabel('Time (min)')
plt.ylabel('Output')
plt.title('Transfer Function Response (step at 150 min)')
plt.grid(which='both', linestyle='--', alpha=0.7)

# Highlight controller activation
plt.axvline(150, color='r', linestyle='--', label='Controller ON (150 min)')

# Adjust x-axis ticks and limits for better visibility
plt.xticks(np.arange(140, 171, 2))   # focus around 150 min
plt.xlim(148, 160)

plt.legend()
plt.minorticks_on()
# plt.show()



##### find the inverse of the transfer functions:
import sympy as sp

def inverse_laplace_transform(eq, s):
    t = sp.symbols('t', real=True)
    output = sp.inverse_laplace_transform(eq, s, t)
    return output  # return it instead of print

# define symbols
lam, s, t = sp.symbols('lam s t', real=True, positive=True)

# define the Laplace-domain expression (note parentheses)
equation_term1 = (lam**2 * s) / (2*lam*s - 2*lam + 1)
equation_term2 = (-2*lam+1)/(lam**2*s**2-2*lam*s+1)
# compute inverse Laplace in terms of lam
inv_eq = inverse_laplace_transform(equation_term1, s)
inv_eq2 = inverse_laplace_transform(equation_term2,s)

print(inv_eq2)




