import ConfigParser


def ReadConfigFile(config_file):
        '''
        Read all the sections in the configuration file.
        @param config_file Configuration file to be readed
        '''
        config_parser = ConfigParser.ConfigParser()
        config_parser.read(config_file)
        
        # Read every section in the file
        sections = config_parser.sections()
        config_options = dict()
        for sect in sections:
            config_options[sect] = ConfigSectionMap(config_parser = config_parser, section = sect)
        
        return config_options

def WriteConfigFile(config_dict, file_name):
        '''
        Write the simulation configuration into a file
        '''
        parser = ConfigParser.ConfigParser()
        
        # Add every section to the file
        for sec in config_dict.keys():
            parser.add_section(sec)
            
            for key in config_dict[sec].keys():
                parser.set(sec, key, config_dict[sec][key])

        with open(file_name, 'w') as f:
            parser.write(f)
        
def ConfigSectionMap(config_parser, section):
    '''
    This function extracts all the properties in a specific section of the file.
    '''
    dict1 = {"name": section}
    options = config_parser.options(section)
    for option in options:
        try:
            dict1[option] = parseOption(config_parser.get(section, option))
        except:
            print 'exception on ',option
            dict1[option] = None
    return dict1
            
def boolify(s):
    if s == 'True' or s == 'true':
            return True
    if s == 'False' or s == 'false':
            return False
    raise ValueError('Not Boolean Value!')

def parseOption(value):
    var = value.split(',')
    new_var = [estimateType(val) for val in var]
        
    if len(new_var)==0:
        return None
    elif len(new_var)==1:
        return new_var[0]
    else:
        return new_var

def estimateType(var):
    '''guesses the str representation of the variables type'''
    var = str(var) #important if the parameters aren't strings...
    for caster in (boolify, int, float):
        try:
            return caster(var)
        except ValueError:
            pass
    
    return var