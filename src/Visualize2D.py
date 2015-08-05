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
    config_dict ['GoC-GoC'] = {'file_name'   : './results/search_state_gocfixed_1p_4p_4c_log.txt',
                              'labels'      : ['GoC-GoC Weight (log10)'],
                              'x'           : 'numpy.log10(numpy.array(abs(param1_values)))'}
        
    # Selected config
    config = 'GoC-GoC'
    
    selected_dict = config_dict[config]
    
    # Load the file with numpy
    data = numpy.loadtxt(selected_dict['file_name'])
    
    # Order the data according to column 0
    data = data[numpy.argsort(data[:,0])]
    
    # Generate the labels for each axis    
    labels = selected_dict['labels']    
    
    # Extract data
    param1_values = data[:,0]
    av_MI_values = data[:,1]
    std_MI_values = data[:,2]
    
    # Get maximum indexes
    max_index = numpy.argmax(av_MI_values)
    
    print 'Best configuration found with values (', param1_values[max_index], ') -> MI:', av_MI_values[max_index]
    
    # Create the figure    
    fig = plt.figure()
    #ax = fig.gca(projection='3d')
    ax = fig.gca()
    x = eval(selected_dict['x'])
    #print min(x), max(x)
    #x = param1_values*1e9
    #y = numpy.log(numpy.array(axis[1]))/numpy.log(2)
    #print min(y), max(y)
    z = av_MI_values
    
    #ax.errorbar(x, av_MI_values, yerr=std_MI_values)
    ax.plot(x,av_MI_values,color='#000000')
    ax.fill_between(x, av_MI_values-std_MI_values, av_MI_values+std_MI_values, alpha=0.5, edgecolor='#3F7F4C', facecolor='#7EFF99',linewidth=0.2)
    # Interpolate the data to generate the mesh
    ax.set_title('Mutual Information')
    ax.set_xlabel(labels[0])
    matplotlib.pylab.savefig('mutual_information.svg', format='svg', dpi=1200) 
    matplotlib.pylab.show()
    
    pass

