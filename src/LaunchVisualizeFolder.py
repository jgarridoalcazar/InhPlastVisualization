#! /usr/bin/env python

import glob
import sys
import os.path
import numpy

import scipy.interpolate
import matplotlib.pyplot as plt

def get_parameter_values(orig_name_list):
    '''
    Search the parameter values in the folder names. They have to be specified as blablabla_param1value_param2value
    '''
    param1_values = []
    param2_values = []
    seed_values = []
    name_list = []
    
    for name in orig_name_list:
        split = name.split('_')
        
        if len(split)<4:
            print 'Warning:',name,'can not be splitted. Pass'
        else:
            seed_values.append(float(split[-3]))
            param1_values.append(float(split[-2]))
            param2_values.append(float(split[-1]))
            name_list.append(name)
    
    return name_list, seed_values, param1_values, param2_values    
    
        
    
def search_folder_names(name_pattern):
    '''
    Search all the folder names that match the specified pattern.
    '''
    return [name for name in glob.glob(name_pattern) if os.path.isdir(name)]
    
    
def visualize_results(name_list, seed_values, param1_values, param2_values):
    '''
    Visualize the results only if the number of parameters is 1 or 2. For the moment it only works
    with two parameters, one pattern and one cell.
    '''

    
    # Generate the labels for each axis    
    labels = ['Param1','Param2']    
    axis = [param1_values,param2_values]
    
    best_MI = 0.0
    best_config = (0,0)
    
    MI_Values = dict()
    Freq_Values = dict()
    
    for param1 in param1_values:
        for param2 in param2_values:
            MI_Values[param1,param2] = []
            Freq_Values[param1,param2] = []
    
    # Extract every mutual information to explore
    for index, name in enumerate(name_list):        
        # Load the MI file
        try:
            mi_file = name + '/' + 'mutual_information'
            value = numpy.loadtxt(mi_file)
            print 'Loaded value in file', mi_file,':', value
            
            MI_Values[param1_values[index],param2_values[index]].append(value[0])
            Freq_Values[param1_values[index],param2_values[index]].append(value[1])
            
        except IOError:
            print "Warning:",mi_file,'does not exist. Using default value 0'
            
            MI_Values[param1_values[index],param2_values[index]].append(0.)
            Freq_Values[param1_values[index],param2_values[index]].append(0.)
            
    MI_Av_Values = []
    Freq_Av_Values =[]
    
    for param1, param2 in zip(param1_values,param2_values):
        Average_MI = numpy.average(MI_Values[param1,param2])
        MI_Av_Values.append(Average_MI)
        
        Average_Freq = numpy.average(Freq_Values[param1,param2])
        Freq_Av_Values.append(Average_Freq)
        
        
        if Average_MI>best_MI:
            best_MI = Average_MI
            best_config = (param1, param2)
    
    print 'Best configuration found with values', best_config, ' -> MI:', best_MI
    
    # Create the figure    
    fig = plt.figure()
    #ax = fig.gca(projection='3d')
    ax = fig.gca()
    x = numpy.array(axis[0])*1.e9
    y = numpy.array(axis[1])*1.
    z = MI_Av_Values
    
    # Set up a regular grid of interpolation points
    xi, yi = numpy.linspace(min(x), max(x), 100), numpy.linspace(min(y), max(y), 100)
    xi, yi = numpy.meshgrid(xi, yi)

    # Interpolate; there's also method='cubic' for 2-D data such as here
    zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method='linear')

    surf = ax.imshow(zi, vmin=min(z), vmax=max(z), origin='lower', extent=[min(x), max(x), min(y), max(y)])
    #surf = ax.plot_surface(xi, yi, zi, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.scatter(x, y, c=z)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    # Interpolate the data to generate the mesh
    ax.set_title('Mutual Information')
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    fig.savefig('mutual_information.png')
    
    # Create the figure    
    fig = plt.figure()
    #ax = fig.gca(projection='3d')
    ax = fig.gca()
    x = numpy.array(axis[0])*1.e9
    y = numpy.array(axis[1])*1.
    z = Freq_Av_Values
    
    # Set up a regular grid of interpolation points
    xi, yi = numpy.linspace(min(x), max(x), 100), numpy.linspace(min(y), max(y), 100)
    xi, yi = numpy.meshgrid(xi, yi)

    # Interpolate; there's also method='cubic' for 2-D data such as here
    zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method='linear')

    surf = ax.imshow(zi, vmin=min(z), vmax=max(z), origin='lower', extent=[min(x), max(x), min(y), max(y)])
    #surf = ax.plot_surface(xi, yi, zi, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.scatter(x, y, c=z)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    # Interpolate the data to generate the mesh
    ax.set_title('Average Firing Frequency')
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    fig.savefig('firing_frequency.png')
        
    #plt.show()

if __name__ == "__main__":
    
    if len(sys.argv)==1:
        print sys.argv
        print 'Error: File name pattern has not been specified. Usage:',sys.argv[0],'file_name_pattern'
        sys.exit(1)
    
    file_name_pattern = sys.argv[1]
    
    file_names = search_folder_names(file_name_pattern)
    
    selected_names, seed_values, param1_values, param2_values = get_parameter_values(file_names)
    
    visualize_results(selected_names, seed_values, param1_values, param2_values)
    
    pass

