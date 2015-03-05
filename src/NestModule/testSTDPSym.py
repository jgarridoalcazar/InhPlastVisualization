import nest
import numpy as np
import pylab as pl
import math

try:
    nest.Install('glplasticitymodule')
except nest.NESTError:
    print('NEST Error caught on loading user module. Retrying...')
    nest.Install('glplasticitymodule')

nest.set_verbosity('M_ERROR')
        

nest.ResetKernel()
nest.SetKernelStatus({"resolution": 0.1, "print_time": False})

# display recordables for illustration
# print('iaf_cond_exp_ip recordables: {0}'.format(nest.GetDefaults('iaf_cond_exp_ip')['recordables']))
print('iaf_cond_exp_sym recordables: {0}'.format(nest.GetDefaults('iaf_cond_exp_sym')['recordables']))

# create neuron and multimeternest.Connect(g_variable, n, sym_spec=syn_dict)nest.Connect(g_variable, n, sym_spec=syn_dict)
n_spikes = 5
tau_sym = 20.0
first_spike = 100.0
spike_offset = 5.0
inter_spike = 50.0

simulation_time = 0.5e3
weight_step = 1

# n = nest.Create('iaf_cond_exp_ip', 2, params = {'tau_syn_ex': 1.0, 
#                                              'V_reset': -65.0, 
#                                              'E_L': -65.0, 
#                                              'E_ex': 0.0, 
#                                              'E_in': -80.0, 
#                                              'V_th': -50.0,
#                                              'beta':1.2, 
#                                              'tau_ip': 1.0e99, 
#                                              'epsilon_rC': 42.0, 
#                                              'epsilon_rR':600.0,
#                                              'tau_sym': tau_sym})
# nest.SetStatus(n,{'r_C':1.0,'g_L':0.6, 'I_e':0.0})
n = nest.Create('iaf_cond_exp_sym', 2, params = {'tau_syn_ex': 1.0, 
                                             'V_reset': -65.0, 
                                             'E_L': -65.0, 
                                             'E_ex': 0.0, 
                                             'E_in': -80.0, 
                                             'V_th': -50.0,
                                             'C_m': 250.0,
                                             'g_L': 16.667,
                                             'tau_sym': tau_sym})
nest.SetStatus(n,{'V_m':-65.0})

# n = nest.Create('iaf_cond_exp', 2)

#g1 = nest.Create('sinusoidal_poisson_generator', n=1, params={'dc': 200.0, 'ac': 200.0,'freq': 0.5, 'phi': 0.0})
#g2 = nest.Create('sinusoidal_poisson_generator', n=1, params={'dc': 200.0, 'ac': 200.0, 'freq': 0.5, 'phi': math.pi})
g1spikes = [first_spike + index*inter_spike for index in range(0,n_spikes)]
g2spikes = [first_spike + spike_offset + index*inter_spike for index in range(0,n_spikes)]
g1spikes.append(simulation_time-5)
g2spikes.append(simulation_time-5)
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
nest.Connect([n[1]], [n[0]], syn_spec='my_stdp')
            
nest.CopyModel('static_synapse','excitatory',
            {'weight':100.0, 
            'delay':1.0})
nest.Connect(g1, [n[0]], syn_spec='excitatory')
nest.Connect(g2, [n[1]], syn_spec='excitatory')

# Get initial weights
stdp_conn1 = nest.GetConnections(source=[n[0]],target=[n[1]])
stdp_conn2 = nest.GetConnections(source=[n[1]],target=[n[0]])
initial_weights1 = nest.GetStatus(stdp_conn1,'weight')
initial_weights2 = nest.GetStatus(stdp_conn2,'weight')
print 'Initial weights:', initial_weights1, initial_weights2

weight_register1 = [initial_weights1]
weight_register2 = [initial_weights2]

for i in np.arange(0.0,simulation_time,weight_step):
    nest.Simulate(weight_step)
    weight_register1.append(nest.GetStatus(stdp_conn1, 'weight'))
    weight_register2.append(nest.GetStatus(stdp_conn2, 'weight'))
    
# simulate

final_weights1 = nest.GetStatus(stdp_conn1, 'weight')
final_weights2 = nest.GetStatus(stdp_conn2, 'weight')
print 'Final weights:', final_weights1, final_weights2

# obtain and display data
vmtime = nest.GetStatus(m)[0]['events']['times']
vmcell = nest.GetStatus(m)[0]['events']['senders']
vmvalue = nest.GetStatus(m)[0]['events']['V_m']

pl.clf()

sp = nest.GetStatus(spikedetector)[0]['events']['times']
cell = nest.GetStatus(spikedetector)[0]['events']['senders']

time_bin = 2.

pl.subplot(411)
h, e = np.histogram(sp[cell==g1], bins=np.arange(0., simulation_time+1., time_bin))
pl.step(e[:-1], h * 1000 / time_bin, where='post', label='Stim 1')
h, e = np.histogram(sp[cell==g2], bins=np.arange(0., simulation_time+1., time_bin))
pl.step(e[:-1], h * 1000 / time_bin, where='post', label='Stim 2')
pl.ylabel('Spike/s')
pl.legend()

pl.subplot(412)
h, e = np.histogram(sp[cell==n[0]], bins=np.arange(0., simulation_time+1., time_bin))
pl.step(e[:-1], h * 1000 / time_bin, where='post', label='Cell 0')
h, e = np.histogram(sp[cell==n[1]], bins=np.arange(0., simulation_time+1., time_bin))
pl.step(e[:-1], h * 1000 / time_bin, where='post', label='Cell 1')
pl.ylabel('Spike/s')
pl.legend()

pl.subplot(413)
pl.plot(vmtime[vmcell==n[0]], vmvalue[vmcell==n[0]], label='Cell 0')
pl.plot(vmtime[vmcell==n[1]], vmvalue[vmcell==n[1]], label='Cell 1')
pl.ylabel('V_m (mV)')
pl.xlabel('Time [ms]')
pl.legend()

pl.subplot(414)
pl.plot(np.arange(0.0,simulation_time+weight_step,weight_step), weight_register2, label='1-0')
pl.plot(np.arange(0.0,simulation_time+weight_step,weight_step), weight_register1, label='0-1')
pl.ylabel('Weight')
pl.legend()

pl.show()
