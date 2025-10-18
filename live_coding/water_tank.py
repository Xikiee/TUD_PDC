"""coding tutorial 1, water tank

this code will not run as some part are incomplete, however this is a skeleton code of how to model a cstr 
"""

#import
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp



# Parameters 
rho =1000
g = 9.81

h_tank = 2.0 #[m]
F_out = 0.8149 
A_tank = 

k = 
P_tank= 
P_downstream = 
density = 


#Utilities
""" functions, calculation, etc"""

def Fin_profile(t):
    if t<25.0:
        return 0.8145
    elif t<=50.0:
        return 0.4
    else:
        return 0.9

### Model ###

def water_level_ode(t,y):
    h_liq = y[0]

    F_out = k * ((density*g*h_liq + P_tank -P_downstream)**0.5)/density 

    dhdt = (Fin_profile(t) - F_out)/A_tank

    # physical boundary 

    h_liq = np.clip(h_liq,0.0,)
    if h_liq <=0 and dhdt <=0:
        dhdt = 0
    if h_liq >= h_tank and dhdt>=0:
        dhdt = 0

    return np.array(dhdt)



### integration ###

tstart = 0.0
tend = 100.0
trange= (tstart,tend)
teval = np.linspace(tstart,tend,2000)


solution = solve_ivp(
    fun = water_level_ode,
    t_span= trange,
    y0 = y0,
    method = "Radau",
    t_eval = teval,
    rtol = 1e-3,
    atol=1e-6,
)

### Results and intepretation 

plt.figure()
plt.plot(teval, solution.y[0], label = "solution" )
plt.show()
