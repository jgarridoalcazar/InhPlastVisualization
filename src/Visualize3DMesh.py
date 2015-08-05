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
    config_dict ['MF-GoC'] = {'file_name'   : './results/search_state_mfgoc_2p_4p_4c_NoNormIP_40mf.txt',
                              'labels'      : ['MF-GoC Max Weight (x10^-9)','MF-GoC Ratio'],
                              'x'           : 'param1_values*1e9',
                              'y'           : 'param2_values'}
    config_dict ['GoC-GoC'] = {'file_name'   : './results/search_state_gocgoc_2p_4p_4c_2.txt',
                              'labels'      : ['GoC-GoC Max Weight (log10)','GoC-GoC Ratio'],
                              'x'           : 'numpy.log10(numpy.array(abs(param1_values)))',
                              'y'           : 'param2_values'}
    
    
    # Selected config
    config = 'MF-GoC'
    
    selected_dict = config_dict[config]
    
    # Load the file with numpy
    data = numpy.loadtxt(selected_dict['file_name'])
    
    # Generate the labels for each axis    
    labels = selected_dict['labels']    
    
    # Extract data
    param1_values = data[:,0]
    param2_values = data[:,1]
    av_MI_values = data[:,2]
    std_MI_values = data[:,3]
    
    # Get maximum indexes
    max_index = numpy.argmax(av_MI_values)
    
    print 'Best configuration found with values (', param1_values[max_index], ',', param2_values[max_index], ') -> MI:', av_MI_values[max_index]
    
    # Create the figure    
    fig = plt.figure()
    #ax = fig.gca(projection='3d')
    ax = fig.gca()
    x = eval(selected_dict['x'])
    #print min(x), max(x)
    #x = param1_values*1e9
    #y = numpy.log(numpy.array(axis[1]))/numpy.log(2)
    #print min(y), max(y)
    y = eval(selected_dict['y'])
    z = av_MI_values
    
    # Set up a regular grid of interpolation points
    xi, yi = numpy.linspace(min(x), max(x), 100), numpy.linspace(min(y), max(y), 100)
    xi, yi = numpy.meshgrid(xi, yi)
    
    # Interpolate; there's also method='cubic' for 2-D data such as here
    zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method='linear')

    surf = ax.imshow(zi, origin='lower', extent=[min(x), max(x), min(y), max(y)],  vmin=0, vmax=0.6)
    #surf = ax.plot_surface(xi, yi, zi, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.scatter(x, y, c=z, vmin=0, vmax=0.6)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    # Interpolate the data to generate the mesh
    ax.set_title('Mutual Information')
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    matplotlib.pylab.savefig('mutual_information.svg', format='svg', dpi=1200) 
    matplotlib.pylab.show()
    
    # Create the figure    
    fig = plt.figure()
    #ax = fig.gca(projection='3d')
    ax = fig.gca()
    z = std_MI_values
    
    # Interpolate; there's also method='cubic' for 2-D data such as here
    zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method='linear')

    surf = ax.imshow(zi, origin='lower', extent=[min(x), max(x), min(y), max(y)], vmin=0, vmax=0.2)
    #surf = ax.plot_surface(xi, yi, zi, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.scatter(x, y, c=z, vmin=0, vmax=0.2)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    # Interpolate the data to generate the mesh
    ax.set_title('Standard Desviation')
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    matplotlib.pylab.savefig('std_mutual.svg', format='svg', dpi=1200) 
    matplotlib.pylab.show()
    
    # Create the figure    
    fig = plt.figure()
    #ax = fig.gca(projection='3d')
    ax = fig.gca()
    z = numpy.zeros(av_MI_values.shape)
    z[av_MI_values>0.0] = std_MI_values[av_MI_values>0.0]/av_MI_values[av_MI_values>0.0]
    z[av_MI_values==0.0] = max(z)
    
    # Interpolate; there's also method='cubic' for 2-D data such as here
    zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method='linear')

    surf = ax.imshow(zi, origin='lower', extent=[min(x), max(x), min(y), max(y)], vmin=0, vmax=1)
    #surf = ax.plot_surface(xi, yi, zi, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.scatter(x, y, c=z, vmin=0, vmax=1)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    # Interpolate the data to generate the mesh
    ax.set_title('CV')
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    matplotlib.pylab.savefig('cv_mutual.svg', format='svg', dpi=1200) 
    matplotlib.pylab.show()

    
    pass

