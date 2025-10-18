import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import scipy 
import control as ctrl

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



lam = 1
gain1_const = R*80/(M_H2*V1gas) 
gain2_const = R*80/(M_O2*V2gas)
print(f"gain for hydrogen {gain1_const}, gain for oxygen {gain2_const}")