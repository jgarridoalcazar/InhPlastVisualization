import nest
import numpy as np
import pylab as pl

try:
    nest.Install('glplasticitymodule')
except nest.NESTError:
    print('NEST Error caught on loading user module. Retrying...')
    nest.Install('glplasticitymodule')

nest.ResetKernel()
nest.SetKernelStatus({"resolution": 0.1, "print_time": True})


# display recordables for illustration
print('iaf_cond_exp_ip recordables: {0}'.format(nest.GetDefaults('iaf_cond_exp_ip')['recordables']))

# create neuron and multimeter
n = nest.Create('iaf_cond_exp_ip', 
                params = {'tau_syn_ex': 1.0, 'V_reset': -65.0, 'E_L': -65.0, 'E_ex': 0.0, 'E_in': -80.0, 'V_th': -50.0,
                           'beta':0.8, 'tau_ip': 1200e3, 'epsilon_rC': 34.55, 'epsilon_rR':4923.8826})
nest.SetStatus(n,{'r_C':0.3,'g_L':6.0, 'I_e':0.0})

g = nest.Create('sinusoidal_poisson_generator', n=1, params={'dc': 50.0, 'ac': 0.0,
                                                              'freq': 0.1, 'phi': 0.0})

m = nest.Create('multimeter',
                params = {'withtime': True, 
                          'interval': 0.1,
                          'record_from': ['V_m', 'g_L', 'r_C']})

spikedetector = nest.Create("spike_detector")
nest.Connect(n, spikedetector)
nest.Connect(g, spikedetector)

nest.Connect(m, n)
nest.Connect(g, n)
simulation_time = 4.0e6
# simulate
nest.Simulate(simulation_time/2)

# Change poisson generator frequency
nest.SetStatus(g,{'dc':200.0,'ac':10.0})
nest.Simulate(simulation_time/2)

# obtain and display data
events = nest.GetStatus(m)[0]['events']
t = events['times'];

pl.clf()

# pl.subplot(411)
# pl.plot(t, events['V_m'])
# pl.ylabel('Membrane potential [mV]')
# pl.xlabel('Time [ms]')

time_bin = 5000.
spike_activity = nest.GetStatus(spikedetector)[0]['events']
spoutput = np.array([time for time,cellid in zip(spike_activity['times'],spike_activity['senders']) if cellid==n[0]])
spinput = np.array([time for time,cellid in zip(spike_activity['times'],spike_activity['senders']) if cellid==g[0]])

avFiringRateInit1 = np.count_nonzero(spoutput<10000.)/10.
aux1 = np.array(spoutput>(simulation_time/2-10000.))
aux2 = np.array(spoutput<(simulation_time/2))
avFiringRateEnd1 = np.count_nonzero(aux1 & aux2)/10.
aux1 = np.array(spoutput>(simulation_time/2))
aux2 = np.array(spoutput<(simulation_time/2+10000.))
avFiringRateInit2 = np.count_nonzero(aux1 & aux2)/10.
avFiringRateEnd2 = np.count_nonzero(spoutput>(simulation_time-10000.))/10.
print 'Average firing rate in 10 seconds. Input 1: Init -', avFiringRateInit1, '. End - ', avFiringRateEnd1, '. Input 2: Init - ', avFiringRateInit2, '. End -', avFiringRateEnd2

pl.subplot(411)
h, e = np.histogram(spinput, bins=np.arange(0., simulation_time+1., time_bin))
pl.plot(e[:-1]/1000.0,h*1000/time_bin)
pl.title('Input activity')
pl.ylabel('Spikes per second')

pl.subplot(412)
h, e = np.histogram(spoutput, bins=np.arange(0., simulation_time+1., time_bin))
#pl.step(e[:-1], h * 1000 / time_bin, where='post')
pl.plot(e[:-1]/1000.0,h*1000/time_bin)
pl.title('Output activity')
pl.ylabel('Spikes per second')

pl.subplot(413)
pl.plot(t/1000.0, events['g_L'])
pl.xlabel('Time [ms]')
pl.ylabel('Leak conductance [nS]')

pl.subplot(414)
pl.plot(t/1000.0, 1./events['r_C'])
pl.xlabel('Time [ms]')
pl.ylabel('Cm [pF]')

pl.savefig('intrinsicplasticity.svg', format='svg', dpi=1200)
pl.show()
 
