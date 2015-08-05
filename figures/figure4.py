#! /usr/bin/env python

import glob
import sys
import os.path
import numpy

import scipy.interpolate
import matplotlib.pyplot as plt
import matplotlib.pylab
    
if __name__ == "__main__":
    
    config_dict = {}
    config_dict ['iSTDP'] = {'MIRatio' : numpy.array([0.61747, 0.54977, 0.526058, 0.55035]),
                             'StdMI' : numpy.array([0.0883276, 0.13493, 0.0982955, 0.0944069])}
    config_dict ['NoiSTDP'] = {'MIRatio' : numpy.array([0.59872, 0.515379, 0.4934648, 0.515379]),
                             'StdMI' : numpy.array([0.17072, 0.15975, 0.1049667, 0.082697])}
    config_dict ['NoInh'] = {'MIRatio' : numpy.array([0.3287918, 0.383137, 0.398203, 0.382614]),
                             'StdMI' : numpy.array([0.324373, 0.168354, 0.1232293, 0.121337])}
    
    NumCells = [1,2,4,8]
    
        
    # Create the figure    
    fig = plt.figure()
    #ax = fig.gca(projection='3d')
    ax = fig.gca()
    #ax.errorbar(NumCells, config_dict['iSTDP']['MIRatio'], yerr=config_dict['iSTDP']['StdMI'])
    ax.plot(NumCells,config_dict['iSTDP']['MIRatio'],color='#99FF99',marker='o',markersize=10)
    ax.fill_between(NumCells, config_dict['iSTDP']['MIRatio']-config_dict['iSTDP']['StdMI'], config_dict['iSTDP']['MIRatio']+config_dict['iSTDP']['StdMI'], alpha=0.5, edgecolor='#99FF99', facecolor='#99FF99',linewidth=0.2)
    #ax.errorbar(NumCells, config_dict['NoiSTDP']['MIRatio'], yerr=config_dict['NoiSTDP']['StdMI'])
    ax.plot(NumCells,config_dict['NoiSTDP']['MIRatio'],color='#FF9999',marker='^',markersize=10)
    ax.fill_between(NumCells, config_dict['NoiSTDP']['MIRatio']-config_dict['NoiSTDP']['StdMI'], config_dict['NoiSTDP']['MIRatio']+config_dict['NoiSTDP']['StdMI'], alpha=0.5, edgecolor='#FF9999', facecolor='#FF9999',linewidth=0.2)
    #ax.errorbar(NumCells, config_dict['NoInh']['MIRatio'], yerr=config_dict['NoInh']['StdMI'])
    ax.plot(NumCells,config_dict['NoInh']['MIRatio'],color='#9999FF',marker='s',markersize=10)
    ax.fill_between(NumCells, config_dict['NoInh']['MIRatio']-config_dict['NoInh']['StdMI'], config_dict['NoInh']['MIRatio']+config_dict['NoInh']['StdMI'], alpha=0.5, edgecolor='#9999FF', facecolor='#9999FF',linewidth=0.2)
    # Interpolate the data to generate the mesh
    ax.legend(['iSTDP', 'Fixed Inh.', 'NoInh'])
    ax.set_title('Mutual Information')
    ax.set_xlabel('Num. Cells and Patterns')
    matplotlib.pylab.savefig('mutual_information_patterns.svg', format='svg', dpi=1200) 
    matplotlib.pylab.show()
    
    config_dict = {}
    config_dict ['iSTDP'] = {'MIRatio' : numpy.array([0.526058, 0.627795, 0.68175246, 0.721692, 0.7562833, 0.7824016, 0.80805]),
                             'StdMI' : numpy.array([0.0982955, 0.08392936, 0.0694552, 0.0544185, 0.0476758, 0.052868, 0.0514911])}
    config_dict ['NoiSTDP'] = {'MIRatio' : numpy.array([0.4934648, 0.572022, 0.62083862, 0.6452358, 0.678271, 0.7100384, 0.740874313]),
                             'StdMI' : numpy.array([0.1049667, 0.0976374, 0.0924136, 0.085242, 0.07762668, 0.0732983, 0.0763218])}
    config_dict ['NoInh'] = {'MIRatio' : numpy.array([0.3982035, 0.4784349, 0.526588, 0.571099, 0.60192, 0.62437, 0.64089156]),
                             'StdMI' : numpy.array([0.1232293, 0.135739, 0.138676, 0.141548, 0.142392, 0.140747, 0.1381065])}
    
    NumCells = [4,6,8,10,12,14,16]
    
    # Create the figure    
    fig = plt.figure()
    #ax = fig.gca(projection='3d')
    ax = fig.gca()
    #ax.errorbar(NumCells, config_dict['iSTDP']['MIRatio'], yerr=config_dict['iSTDP']['StdMI'])
    ax.plot(NumCells,config_dict['iSTDP']['MIRatio'],color='#99FF99',marker='o',markersize=10)
    ax.fill_between(NumCells, config_dict['iSTDP']['MIRatio']-config_dict['iSTDP']['StdMI'], config_dict['iSTDP']['MIRatio']+config_dict['iSTDP']['StdMI'], alpha=0.5, edgecolor='#99FF99', facecolor='#99FF99',linewidth=0.2)
    #ax.errorbar(NumCells, config_dict['NoiSTDP']['MIRatio'], yerr=config_dict['NoiSTDP']['StdMI'])
    ax.plot(NumCells,config_dict['NoiSTDP']['MIRatio'],color='#FF9999',marker='^',markersize=10)
    ax.fill_between(NumCells, config_dict['NoiSTDP']['MIRatio']-config_dict['NoiSTDP']['StdMI'], config_dict['NoiSTDP']['MIRatio']+config_dict['NoiSTDP']['StdMI'], alpha=0.5, edgecolor='#FF9999', facecolor='#FF9999',linewidth=0.2)
    ax.plot(NumCells,config_dict['NoInh']['MIRatio'],color='#9999FF',marker='s',markersize=10)
    ax.fill_between(NumCells, config_dict['NoInh']['MIRatio']-config_dict['NoInh']['StdMI'], config_dict['NoInh']['MIRatio']+config_dict['NoInh']['StdMI'], alpha=0.5, edgecolor='#9999FF', facecolor='#9999FF',linewidth=0.2)
    ax.legend(['iSTDP', 'Fixed Inh.', 'No Inh.'])
    ax.set_title('Mutual Information')
    ax.set_xlabel('Num. Cells and Patterns')
    matplotlib.pylab.savefig('mutual_information_cells.svg', format='svg', dpi=1200) 
    matplotlib.pylab.show()
    
    pass

