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
    config_dict ['Norm10'] = {'MIRatio' : numpy.array([0.0782, 0.4246, 0.5261, 0.5229, 0.4280, 0.4038, 0.4011, 0.4194, 0.3814, 0.3543, 0.3307, 0.3231]),
                             'StdMI' : numpy.array([0.0312, 0.1014, 0.0983, 0.1003, 0.0946, 0.0954, 0.0899, 0.0897, 0.0768, 0.0692, 0.0713, 0.0727])}
    config_dict ['NoNorm10'] = {'MIRatio' : numpy.array([0.0587, 0.2494, 0.5046, 0.2980, 0.2619, 0.2467, 0.2012, 0.1838, 0.2170, 0.2723, 0.3078, 0.3193]),
                             'StdMI' : numpy.array([0.0212, 0.0620, 0.1069, 0.0711, 0.0708, 0.0888, 0.1061, 0.1074, 0.1067, 0.0839, 0.0662, 0.0606])}
    config_dict ['Norm40'] = {'MIRatio' : numpy.array([0.0382, 0.1250, 0.3199, 0.5016, 0.5386, 0.5091, 0.4600, 0.4387, 0.4156, 0.3814, 0.3508, 0.3386]),
                             'StdMI' : numpy.array([0.0042, 0.0446, 0.0549, 0.0908, 0.1020, 0.1254, 0.1268, 0.1166, 0.1039, 0.0822, 0.0740, 0.0743])}
    config_dict ['NoNorm40'] = {'MIRatio' : numpy.array([0.0677, 0.0808, 0.2186, 0.4333, 0.5325, 0.5418, 0.5326, 0.4232, 0.3458, 0.3273, 0.3274, 0.3342]),
                             'StdMI' : numpy.array([0.0080, 0.0232, 0.0555, 0.0794, 0.0825, 0.0899, 0.1114, 0.0922, 0.0803, 0.0863, 0.0920, 0.0910])}
            
    NumCells = [1,5,10,20,30,40,50,60,70,80,90,100]
    
        
    # Create the figure    
    fig = plt.figure()
    #ax = fig.gca(projection='3d')
    ax = fig.gca()
    #ax.errorbar(NumCells, config_dict['iSTDP']['MIRatio'], yerr=config_dict['iSTDP']['StdMI'])
    ax.plot(NumCells,config_dict['Norm10']['MIRatio'],color='#99FF99',marker='o',markersize=10)
    ax.fill_between(NumCells, config_dict['Norm10']['MIRatio']-config_dict['Norm10']['StdMI'], config_dict['Norm10']['MIRatio']+config_dict['Norm10']['StdMI'], alpha=0.5, edgecolor='#99FF99', facecolor='#99FF99',linewidth=0.2)
    #ax.errorbar(NumCells, config_dict['NoiSTDP']['MIRatio'], yerr=config_dict['NoiSTDP']['StdMI'])
    ax.plot(NumCells,config_dict['NoNorm10']['MIRatio'],color='#FF9999',marker='^',markersize=10)
    ax.fill_between(NumCells, config_dict['NoNorm10']['MIRatio']-config_dict['NoNorm10']['StdMI'], config_dict['NoNorm10']['MIRatio']+config_dict['NoNorm10']['StdMI'], alpha=0.5, edgecolor='#FF9999', facecolor='#FF9999',linewidth=0.2)
    # Interpolate the data to generate the mesh
    ax.legend(['Control', 'No Norm. No IP'])
    ax.set_title('Mutual Information')
    ax.set_xlabel('Num. Cells and Patterns')
    matplotlib.pylab.savefig('mutual_information_norm_10.svg', format='svg', dpi=1200) 
    matplotlib.pylab.show()
        
        
    # Create the figure    
    fig = plt.figure()
    #ax = fig.gca(projection='3d')
    ax = fig.gca()
    #ax.errorbar(NumCells, config_dict['iSTDP']['MIRatio'], yerr=config_dict['iSTDP']['StdMI'])
    ax.plot(NumCells,config_dict['Norm40']['MIRatio'],color='#99FF99',marker='o',markersize=10)
    ax.fill_between(NumCells, config_dict['Norm40']['MIRatio']-config_dict['Norm40']['StdMI'], config_dict['Norm40']['MIRatio']+config_dict['Norm40']['StdMI'], alpha=0.5, edgecolor='#99FF99', facecolor='#99FF99',linewidth=0.2)
    #ax.errorbar(NumCells, config_dict['NoiSTDP']['MIRatio'], yerr=config_dict['NoiSTDP']['StdMI'])
    ax.plot(NumCells,config_dict['NoNorm40']['MIRatio'],color='#FF9999',marker='^',markersize=10)
    ax.fill_between(NumCells, config_dict['NoNorm40']['MIRatio']-config_dict['NoNorm40']['StdMI'], config_dict['NoNorm40']['MIRatio']+config_dict['NoNorm40']['StdMI'], alpha=0.5, edgecolor='#FF9999', facecolor='#FF9999',linewidth=0.2)
    # Interpolate the data to generate the mesh
    ax.legend(['Control', 'No Norm. No IP'])
    ax.set_title('Mutual Information')
    ax.set_xlabel('Num. Cells and Patterns')
    matplotlib.pylab.savefig('mutual_information_norm_40.svg', format='svg', dpi=1200) 
    matplotlib.pylab.show()
    
    pass

