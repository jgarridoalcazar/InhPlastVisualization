#include "uego.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fitness.c"
#include <vector>


#include "uego.h"
#include "time.h"

////////////////////////////////////////////////////////////
// $Id: nrealval.cc,v 2.5 1998/03/17 23:14:53 jelasity Exp $
// nrealval.cc
// definition of NDimRealElement::Value()
// objective function library for n dimensional real spaces
////////////////////////////////////////////////////////////
// modification history:
//	Pilar M. Ortigosa 98 01 27 lot of new functions
//	Pilar M. Ortigosa 98 02 10 minor modifications
//	Pilar M. Ortigosa 02 06 09 delays in test functions
////////////////////////////////////////////////////////////



//------------------------------------------------------------------


/*
double	NDimRealElement::Value() {
// INI.Fnum() can be used to select a funcion
// INI.Lowb(i) and INI.Upb(i) are lower and upper bounds on x[i]
// INI.Param(i) is the ith constant parameter, INI.ParamNum() is the
// number of available parameters
// dim is the dimension of space

// --- Variables para el retardo
 clock_t 		rel1,rel2;
 double 	temp;
 int 		count;
 double 	res;
 double 	y,z;
 //---------------------------



	// --- check for bounds ------------------------
	for(long j = 0; j < dim; ++j )
	{
		if( x[j] < INI.Lowb(j) )
		{
			message((char*)"Variable out of range, correcting.",MSG_INFORMATION);
			x[j] = INI.Lowb(j);
		};
		if( x[j] > INI.Upb(j) )
		{
			message((char*)"Variable out of range, correcting.",MSG_INFORMATION);
			x[j] = INI.Upb(j);
		};
	};
	
	
	//-- Los retardos----------------------
	rel1=clock();
	rel2=clock();
 	res=(double)(rel2-rel1)/(double)CLOCKS_PER_SEC;
	 while (res< 0.003) {
		rel2=clock();
		res=(double)(rel2-rel1)/(double)CLOCKS_PER_SEC;
 	};


	 // Levy 3 [-10,10]^2
	//N∫ 16 tesis
	y=0.0;
	for( long j=1; j<6; ++j ) y += j * cos((j-1)*x[0]+j);
	z=0.0;
	for( long j=1; j<6; ++j ) z += j * cos((j+1)*x[1]+j);
	return -y*z;


	//Beale

	//x[0] = 3;
	//x[1] = 0.5;
	//f = 0;

	y = 0.0;
	y = pow((1.5 -x[0]+ x[0]*x[1]),2.0)+pow((2.25-x[0]+x[0]*x[1]*x[1]),2.0)+pow((2.625-x[0]+x[0]*x[1]*x[1]*x[1]),2.0);
	//printf(" %lf", y);
	//getchar();
	return -y;

};
*/

// INI.Lowb(i) and INI.Upb(i) are lower and upper bounds on x[i] max(i)=dim-1;
double	NDimRealElement::Value()
{

std::vector<std::string> param_names(INI.Dimension());
std::vector<double> param_values(INI.Dimension());
double  med_fobj, sigma,varObj;
int count = 1;
double *objFunc=new double[count];
double finalFitness;
unsigned int seed = 12345;

//---Allocation memory
for(int i=0;i<INI.Dimension();i++){
	param_names[i]=INI.ParameterName(i);
	param_values[i]=x[i];
}


//---Computing objective function value
 med_fobj = 0;
for(int i=0;i<count;i++){
	objFunc[i] = fitness(seed,param_names, param_values);
	//objFunc[i] = i;
	med_fobj += objFunc[i];
	seed++;
};

//---Computing average fitness value
 med_fobj /=count;

//---Computing sigma
varObj=0;
for (long i=0;i<count;i++)
 	varObj += pow((double)objFunc[i] - med_fobj,2.0);
 	
varObj /= (double)count;
sigma = sqrt(varObj);

//finalFitness = med_fobj/sigma;
finalFitness = med_fobj;
return finalFitness;
}

