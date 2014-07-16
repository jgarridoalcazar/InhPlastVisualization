import nest
import numpy as np
import pylab as pl

nest.Install('glplasticitymodule')

nest.ResetKernel()
nest.SetKernelStatus({"resolution": 0.1, "print_time": True})


# display recordables for illustration
print('iaf_cond_exp_ip recordables: {0}'.format(nest.GetDefaults('iaf_cond_exp_ip')['recordables']))

# create neuron and multimeter
n = nest.Create('iaf_cond_exp_ip', 
                params = {'tau_syn_ex': 1.0, 'V_reset': -65.0, 'E_L': -65.0, 'E_ex': 0.0, 'E_in': -80.0, 'V_th': -50.0,
                           'beta':1.2, 'tau_ip': 1.0e6, 'epsilon_rC': 42.0, 'epsilon_rR':600.0})
nest.SetStatus(n,{'r_C':1.0,'g_L':0.6, 'I_e':0.0})

g = nest.Create('sinusoidal_poisson_generator', n=1, params={'dc': 100.0, 'ac': 100.0,
                                                              'freq': 0.5, 'phi': 0.0})

m = nest.Create('multimeter',
                params = {'withtime': True, 
                          'interval': 0.1,
                          'record_from': ['V_m', 'g_L', 'r_C']})

spikedetector = nest.Create("spike_detector")
nest.Connect(n, spikedetector)


nest.Connect(m, n)
nest.Connect(g, n)
simulation_time = 1.e6
# simulate
nest.Simulate(simulation_time)

# obtain and display data
events = nest.GetStatus(m)[0]['events']
t = events['times'];

pl.clf()

pl.subplot(411)
pl.plot(t, events['V_m'])
pl.ylabel('Membrane potential [mV]')
pl.xlabel('Time [ms]')

time_bin = 1000.
sp = nest.GetStatus(spikedetector)[0]['events']['times']

avFiringRate = np.count_nonzero(sp>(simulation_time-100000.))/100.
print 'Average firing rate in the last 100seconds:', avFiringRate

pl.subplot(412)
h, e = np.histogram(sp, bins=np.arange(0., simulation_time+1., time_bin))
pl.step(e[:-1], h * 1000 / time_bin, where='post')
pl.title('PST histogram and firing rates')
pl.ylabel('Spikes per second')

pl.subplot(413)
pl.plot(t, events['g_L'])
pl.xlabel('Time [ms]')
pl.ylabel('Leak conductance [nS]')

pl.subplot(414)
pl.plot(t, events['r_C'])
pl.xlabel('Time [ms]')
pl.ylabel('Inverse of Cm [pF-1]')

pl.show()
