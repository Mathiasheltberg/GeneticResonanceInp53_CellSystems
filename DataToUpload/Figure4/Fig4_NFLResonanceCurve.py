from scipy.integrate import RK45,solve_ivp
import numpy as np
from numba import njit, prange
import time
import matplotlib.pyplot as plt
from Functions_Simulations import *

np.random.seed(42)
Amps = []
OM = []
Pars = []
Tmax = 100
ts = 0.01
for test2 in range(1000):
    print(test2)
    A1 = 0.05
    omega = 1.0
    lamR = -.1
    acc = 0
    while acc == 0:
        b2 = A1 +3*np.random.random()
        b3 = 5*np.random.random()
        b4 = 20*np.random.random()
    
        beta = 2*lamR + b3 + 2*b2**2*b4/b3**2
        gamma = 2*lamR*b2**2*b4/b3**2 - b2**2*b4/b3 + b4**2*b2**4/b3**4
        Det = beta**2-4*gamma
        if (Det>0):
            b1 = 0.5*(-beta + np.sqrt(Det))
            if (b1>0):
                acc = 1
                Delta = b3*b1 + b4*b2**2/b3
                om0 = 0.5*np.sqrt(4*Delta - (2*lamR)**2)
    ts = np.linspace(0, Tmax, Tmax*100)
    bru1 =  lambda T,Y: [-b3*Y[0] + b1*Y[1] + b4*Y[0]**2*Y[1],
                          A1*np.sin(omega*T) + b2 - b1*Y[1] - b4*Y[0]**2*Y[1]]
     
    x0 = b2/b3
    y0 = b2*b3**2/(b1*b3**2 + b2**2*b4)
    sol = solve_ivp (bru1, [0, Tmax], [x0+1, y0+1],t_eval=ts,rtol=1e-5)
     
    T = sol.t
    Y = sol.y
    t = T;
    x = Y[:][0];
    y = Y[:][1];

    n_sig= 2
    peaks,valleys = Find_Peaks(x,n_sig)
    Amp = 0
    if (len(peaks)>3 and len(valleys)>3):
        Top = x[peaks[-3:]]
        Bot = x[valleys[-3:]]
        Amp = np.mean(0.5*(Top-Bot))/np.mean(x)

    
    Amps.append(Amp)
    OM.append(om0)

OM = np.array(OM)
Amps = np.array(Amps)

plt.plot(OM,Amps,'xm')
plt.show()
plt.hist(Amps,30,facecolor = 'm',edgecolor = 'k',lw=3)
plt.show()
