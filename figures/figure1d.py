import numpy as np
import pylab as pl
from math import atan, exp, sin, pi

# Excitatory STDP figure

time = np.arange(-0.1,0.1,1e-4)

tau_plus = 16.8e-3
tau_minus = 33.7e-3
minus_plus_ratio = 1.30

estdp = np.zeros(len(time))
estdp[time<0] = -np.exp(time[time<0]/tau_minus)
estdp[time>0] = np.exp(-time[time>0]/tau_plus)/minus_plus_ratio

tau_plus = 125.0e-3
tau_minus= 195.6e-3
Altp = 1
Altd = 0.5*Altp
C = 1./(exp(-atan(pi/2.0)*4./pi)*pow(sin(atan(pi/2)),2))
tau_minus = tau_plus/(atan(pi/2)*2./pi)

time2 = np.arange(-0.3,0.3,1e-4)

istdp = Altp*np.exp(-np.abs(time2)/tau_plus)*np.cos(time2*pi/(tau_plus*2.))**2-Altd*C*np.exp(-2*np.abs(time2)/tau_minus)*np.sin(time2*pi/(tau_minus*2))**2

pl.clf()

pl.subplot(211)
pl.plot(time, estdp)
pl.ylabel('Weight Change (%)')
pl.xlabel('Time difference (s)')

pl.subplot(212)
pl.plot(time2, istdp)
pl.ylabel('Weight Change (%)')
pl.xlabel('Time difference (s)')

pl.savefig('stdp.svg', format='svg', dpi=1200)
pl.show()
 
