import nest
import numpy as np
import pylab as pl

nest.Install('glplasticitymodule')

nest.ResetKernel()
nest.SetKernelStatus({"resolution": 0.1, "print_time": True})


# display recordables for illustration
print('iaf_cond_exp_sto_ip recordables: {0}'.format(nest.GetDefaults('iaf_cond_exp_sto_ip')['recordables']))

# create neuron and multimeter
n = nest.Create('iaf_cond_exp_sto_ip', 
                params = {'t_ref_abs': 2.0, 't_ref': 10.0, 'V_reset': -65.0, 'g_L': 1.0, 'C_m': 50.0, 'E_L': -65.0, 'E_ex': 0.0, 'E_in': -80.0, 'tau_syn_ex': 1.0,
                          'I_e': 0.0, 'ip_rate': 1e-7, 'target_firing': 3.0})
nest.SetStatus(n,{'V_m':-65.0, 'V_th':-65.0, 'r_0':3.0, 'u_alpha': 2.0})

g = nest.Create('sinusoidal_poisson_generator', n=1, params={'dc': 100.0, 'ac': 100.0,
                                                              'freq': 8.0, 'phi': 0.0})

m = nest.Create('multimeter',
                params = {'withtime': True, 
                          'interval': 0.1,
                          'record_from': ['V_m', 'g_ex', 'g_in', 'V_th', 'r_0', 'u_alpha', 'refractoriness', 'gain', 'firing_probability']})

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

pl.subplot(5,2,1)
pl.plot(t, events['V_m'])
pl.ylabel('Membrane potential [mV]')
pl.xlabel('Time [ms]')

pl.subplot(522)
pl.plot(t, events['g_ex'])
pl.xlabel('Time [ms]')
pl.ylabel('Excitatory conductance [nS]')

pl.subplot(523)
pl.plot(t, events['g_in'])
pl.xlabel('Time [ms]')
pl.ylabel('Inhibitory conductance [nS]')

pl.subplot(524)
pl.plot(t, events['V_th'])
pl.xlabel('Time [ms]')
pl.ylabel('Threshold potential [mV]')

pl.subplot(525)
pl.plot(t, events['r_0'])
pl.xlabel('Time [ms]')
pl.ylabel('Firing rate gain [Hz]')

pl.subplot(526)
pl.plot(t, events['u_alpha'])
pl.xlabel('Time [ms]')
pl.ylabel('Potential Scale [mV]')

pl.subplot(527)
pl.plot(t, events['refractoriness'])
pl.xlabel('Time [ms]')
pl.ylabel('Refractoriness')

pl.subplot(528)
pl.plot(t, events['gain'])
pl.xlabel('Time [ms]')
pl.ylabel('Gain')

pl.subplot(529)
pl.plot(t, events['firing_probability'])
pl.xlabel('Time [ms]')
pl.ylabel('Firing Probability')

time_bin = 1000.
sp = nest.GetStatus(spikedetector)[0]['events']['times']
 
avFiringRate = np.count_nonzero(sp>(simulation_time-100000.))/100.
print 'Average firing rate in the last 100seconds:', avFiringRate
 
pl.subplot(5,2,10)
h, e = np.histogram(sp, bins=np.arange(0., simulation_time+1., time_bin))
pl.step(e[:-1], h * 1000 / time_bin, where='post')
pl.title('PST histogram and firing rates')
pl.ylabel('Spikes per second')


pl.show()
