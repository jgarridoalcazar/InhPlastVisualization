from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.pylab
import numpy


import re

file1 = ['./plots/four_patterns/weight-scaling/first-round/test/test-all_4p4c.log',
         './plots/four_patterns/weight-scaling/first-round/test/test-all_4p8c.log',
         './plots/four_patterns/weight-scaling/first-round/test/test-all_4p16c.log']
file2 = ['./plots/four_patterns/weight-scaling/first-round/test/test-NoISTDP_4p4c.log',
         './plots/four_patterns/weight-scaling/first-round/test/test-NoISTDP_4p8c.log',
         './plots/four_patterns/weight-scaling/first-round/test/test-NoISTDP_4p16c.log']
file3 = ['./plots/four_patterns/weight-scaling/first-round/test/test-NoInh_4p4c.log',
         './plots/four_patterns/weight-scaling/first-round/test/test-NoInh_4p8c.log',
         './plots/four_patterns/weight-scaling/first-round/test/test-NoInh_4p16c.log']

def extract_mi_ratio(file):

    print 'Extracting information from',file
    
    output_list = []
    
    
    pattern = r"INFO: Mutual information with seed (\d+) and parameters \[([-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?)\]: ([-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?)"
    
    with open(file) as input:
        for line in input:
            #print 'Reading line:',line
            
            m = re.search(pattern, line)
            if m:
                #print m.group(1),m.group(2),m.group(3)
                try:
                    output_list.append(float(m.group(3)))
                except ValueError:
                    print m.group(3), "is not float"
    
    return output_list

for ind in range(3):
    rvs1 = extract_mi_ratio(file1[ind])
    rvs2 = extract_mi_ratio(file2[ind])
    rvs3 = extract_mi_ratio(file3[ind])

    ttest12, pvalue12 = stats.ttest_ind(rvs1, rvs2, equal_var=False)
    ttest13, pvalue13 = stats.ttest_ind(rvs1, rvs3, equal_var=False)
    ttest23, pvalue23 = stats.ttest_ind(rvs2, rvs3, equal_var=False)
    
    print 'All vs No iSTDP -> ttest:', ttest12, 'pvalue:', pvalue12
    print 'All vs No Inh -> ttest:', ttest13, 'pvalue:', pvalue13
    print 'No iSTDP vs No Inh -> ttest:', ttest23, 'pvalue:', pvalue23
    
    # Create the figure    
    fig = plt.figure()
    #ax = fig.gca(projection='3d')
    plt.subplot(211)
    ax = fig.gca()
        
    #ax.errorbar(x, av_MI_values, yerr=std_MI_values)
    ax.scatter(numpy.ones((1,len(rvs1))),rvs1,color='#00FF00',s=2)
    ax.scatter(numpy.ones((1,len(rvs2)))*2,rvs2,color='#FF0000',s=2)
    ax.scatter(numpy.ones((1,len(rvs3)))*3,rvs3,color='#0000FF',s=2)
    ax.boxplot([rvs1, rvs2, rvs3])
    ax.legend(['iSTDP','Fixed Inh.', 'No Inh.'])
    # Interpolate the data to generate the mesh
    ax.set_title('Mutual Information')
    ax.set_ylabel('Mutual Information Ratio')
    
    hist1, bin_edges1 = numpy.histogram(rvs1, bins=100, range=(0,1), density=False)
    hist2, bin_edges2 = numpy.histogram(rvs2, bins=100, range=(0,1), density=False)
    hist3, bin_edges3 = numpy.histogram(rvs3, bins=100, range=(0,1), density=False)
    
    centers = (bin_edges1[:-1]+bin_edges1[1:])/2
    
    plt.subplot(212)
    ax = fig.gca()
    ax.plot(centers,hist1, color='#00FF00')
    ax.plot(centers,hist2, color='#FF0000')
    ax.plot(centers,hist3, color='#0000FF')
    ax.set_ylabel('Density')
    ax.set_xlabel('Mutual Information Ratio')
    
    matplotlib.pylab.savefig('mutual_information'+str(ind)+'.svg', format='svg', dpi=1200) 
    matplotlib.pylab.show()
    