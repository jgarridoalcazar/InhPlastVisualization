#include "uego.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fitness.c"
#include <vector>

// INI.Lowb(i) and INI.Upb(i) are lower and upper bounds on x[i] max(i)=dim-1;
double	NDimRealElement::Value()
{

double  med_fobj, sigma,varObj;
int count = 5;
double *objFunc=new double[count];
double finalFitness;
std::vector<std::string> param_names(INI.Dimension());
std::vector<double> param_values(INI.Dimension());

//---Setting values
for(int i=0;i<INI.Dimension();i++)  param_names[i] = INI.ParameterName(i);
for(int i=0;i<INI.Dimension();i++)  param_values[i] = x[i];

//---Showing values
/*for(int i=0;i<INI.Dimension();i++)  printf(" %s ", param_names[i].c_str());
for(int i=0;i<INI.Dimension();i++)  printf(" %e ", param_values[i]);
for(int i=0;i<INI.Dimension();i++)  printf(" %e ", x[i]);

getchar();*/

//---Computing objective function value
med_fobj = 0;
for(int i=0;i<count;i++){
	objFunc[i] = fitness("simulation_test",param_names, param_values);
	//objFunc[i] = i;
	med_fobj += objFunc[i];

};

//---Computing average fitness value
 med_fobj /=count;

//---Computing sigma
/*varObj=0;
for (long i=0;i<count;i++)
 	varObj += pow((double)objFunc[i] - med_fobj,2.0);
 	
varObj /= (double)count;
sigma = sqrt(varObj);

finalFitness = med_fobj/sigma;*/
finalFitness= med_fobj;
return finalFitness;

//return(fitness("simulation_test",param_names, param_values));
//getchar();

//delete param_names;
//delete param_values;
}

