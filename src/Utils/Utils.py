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