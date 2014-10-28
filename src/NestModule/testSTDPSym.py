import nest
import numpy as np
import pylab as pl
import math

try:
    nest.Install('glplasticitymodule')
except nest.NESTError:
    print('NEST Error caught on loading user module. Retrying...')
    nest.Install('glplasticitymodule')

nest.sr("M_WARNING setverbosity")
        

nest.ResetKernel()
nest.SetKernelStatus({"resolution": 0.1, "print_time": True})

# display recordables for illustration
print('iaf_cond_exp_ip recordables: {0}'.format(nest.GetDefaults('iaf_cond_exp_ip')['recordables']))

# create neuron and multimeternest.Connect(g_variable, n, sym_spec=syn_dict)nest.Connect(g_variable, n, sym_spec=syn_dict)
n_spikes = 10000
tau_sym = 20.0
first_spike = 100.0
spike_offset = 5.0
inter_spike = 100.0

n = nest.Create('iaf_cond_exp_ip', 2, params = {'tau_syn_ex': 1.0, 
                                             'V_reset': -65.0, 
                                             'E_L': -65.0, 
                                             'E_ex': 0.0, 
                                             'E_in': -80.0, 
                                             'V_th': -50.0,
                                             'beta':1.2, 
                                             'tau_ip': 1.0e99, 
                                             'epsilon_rC': 42.0, 
                                             'epsilon_rR':600.0,
                                             'tau_sym': tau_sym})
nest.SetStatus(n,{'r_C':1.0,'g_L':0.6, 'I_e':0.0})

# n = nest.Create('iaf_cond_exp', 2)

#g1 = nest.Create('sinusoidal_poisson_generator', n=1, params={'dc': 200.0, 'ac': 200.0,'freq': 0.5, 'phi': 0.0})
#g2 = nest.Create('sinusoidal_poisson_generator', n=1, params={'dc': 200.0, 'ac': 200.0, 'freq': 0.5, 'phi': math.pi})
g1spikes = [first_spike + index*inter_spike for index in range(0,n_spikes)]
g2spikes = [first_spike + spike_offset + index*inter_spike for index in range(0,n_spikes)]
g1 = nest.Create('spike_generator', n=1, params={'spike_times': g1spikes})
g2 = nest.Create('spike_generator', n=1, params={'spike_times': g2spikes})

m = nest.Create('multimeter',
                params = {'withtime': True,
                          'interval': 0.1,
                          'record_from': ['V_m']})

spikedetector = nest.Create("spike_detector")
nest.Connect(n, spikedetector, 'all_to_all')
nest.Connect(g1, spikedetector)
nest.Connect(g2, spikedetector)
nest.Connect(m, n, 'all_to_all')

# nest.CopyModel('stdp_synapse_hom','my_stdp',
# 			{'weight':0.0, 
# 			'delay':1.0, 
# 			'Wmax': 10.0, 
# 			'alpha': 0.1, 
# 			'lambda': 0.1})
nest.CopyModel('stdp_sym_synapse_hom','my_stdp',
            {'weight':-0.1, 
            'delay':0.1, 
            'Wmax': -10.0, 
            'alpha': 0.2, 
            'tau_sym':tau_sym, 
            'lambda': 0.01})
           
nest.Connect([n[0]], [n[1]], syn_spec='my_stdp')
            
nest.CopyModel('static_synapse','excitatory',
            {'weight':1.0, 
            'delay':1.0})
nest.Connect(g1, [n[0]], syn_spec='excitatory')
nest.Connect(g2, [n[1]], syn_spec='excitatory')

# Get initial weights
stdp_conn = nest.GetConnections(source=[n[0]],target=[n[1]])
initial_weights = nest.GetStatus(stdp_conn,'weight')
print 'Initial weights:', initial_weights

simulation_time = 3.e4
weight_step = 1

weight_register = [initial_weights]

for i in np.arange(0.0,simulation_time,weight_step):
    nest.Simulate(weight_step)
    weight_register.append(nest.GetStatus(stdp_conn, 'weight'))
    
# simulate

final_weights = nest.GetStatus(stdp_conn, 'weight')
print 'Final weights:', final_weights

# obtain and display data
vmtime = nest.GetStatus(m)[0]['events']['times']
vmcell = nest.GetStatus(m)[0]['events']['senders']
vmvalue = nest.GetStatus(m)[0]['events']['V_m']

pl.clf()

sp = nest.GetStatus(spikedetector)[0]['events']['times']
cell = nest.GetStatus(spikedetector)[0]['events']['senders']

time_bin = 1000.

pl.subplot(421)
h, e = np.histogram(sp[cell==g1], bins=np.arange(0., simulation_time+1., time_bin))
pl.step(e[:-1], h * 1000 / time_bin, where='post')
pl.ylabel('Stim 1 (spike/s)')

pl.subplot(422)
h, e = np.histogram(sp[cell==g2], bins=np.arange(0., simulation_time+1., time_bin))
pl.step(e[:-1], h * 1000 / time_bin, where='post')
pl.ylabel('Stim 2 (spike/s)')

pl.subplot(423)
h, e = np.histogram(sp[cell==n[0]], bins=np.arange(0., simulation_time+1., time_bin))
pl.step(e[:-1], h * 1000 / time_bin, where='post')
pl.ylabel('Cell 0 (spike/s)')

pl.subplot(424)
pl.plot(vmtime[vmcell==n[0]], vmvalue[vmcell==n[0]])
pl.ylabel('V_m cell 0 (mV)')
pl.xlabel('Time [ms]')

pl.subplot(425)
h, e = np.histogram(sp[cell==n[1]], bins=np.arange(0., simulation_time+1., time_bin))
pl.step(e[:-1], h * 1000 / time_bin, where='post')
pl.ylabel('Cell 1 (spike/s)')

pl.subplot(426)
pl.plot(vmtime[vmcell==n[1]], vmvalue[vmcell==n[1]])
pl.ylabel('V_m cell 1 (mV)')
pl.xlabel('Time [ms]')

pl.subplot(427)
pl.plot(np.arange(0.0,simulation_time+weight_step,weight_step), weight_register)
pl.ylabel('Weight')

pl.show()
