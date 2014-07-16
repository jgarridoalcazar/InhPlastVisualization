#! /usr/bin/env python
import glob
import sys
import os
import os.path

def rename_files(folder_list):
    '''
    Search the parameter values in the folder names. They have to be specified as blablabla_param1value_param2value
    '''
    for folder in folder_list:
        file_name_pattern = folder + '/mutual_information_s*'
    
        for name in glob.glob(file_name_pattern):
            if os.path.isfile(name):
                split = name.split('/')
                new_name = '/'.join(split[:-1]) + '/mutual_information'
                print 'Renaming',name,'to',new_name
                os.rename(name, new_name)    
        
    
def search_folder_names(name_pattern):
    '''
    Search all the folder names that match the specified pattern.
    '''
    return [name for name in glob.glob(name_pattern) if os.path.isdir(name)]
    
    

if __name__ == "__main__":
    
    if len(sys.argv)==1:
        print sys.argv
        print 'Error: File name pattern has not been specified. Usage:',sys.argv[0],'file_name_pattern'
        sys.exit(1)
    
    file_name_pattern = sys.argv[1]
    
    file_names = search_folder_names(file_name_pattern)
    
    rename_files(file_names)
    
    pass
