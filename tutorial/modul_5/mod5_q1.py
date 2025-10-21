import numpy as np
import matplotlib.pyplot as plt
import control as ctrl


num = [1]
den = [1,1]
sys = ctrl.TransferFunction(num,den)

# ---- Defining the Delay ----
delay = 0.3
num_delay, den_delay = ctrl.pade(delay, 8)
d_func = ctrl.TransferFunction(num_delay,den_delay)

sys_delay = sys * d_func
time = np.linspace(0,8,8001)
t,x = ctrl.step_response(sys_delay,time)


plt.plot(t,x)
plt.xlabel("x")
plt.ylabel('y')
plt.minorticks_on()
plt.grid(which='both')
plt.xlim(0,8)
plt.show()