#include "uego.h"
#include "time.h"
//#include "internal.h"
#include <stdlib.h>
#include <unistd.h>

//-----------------------------------------------------------------------

#define	MaxLocalOptima 505

double	a[MaxLocalOptima];		// coefficients
	
//------------------------------------------------------------------------

#define opt_(i,j) ( INI.Param( 3 + (i)*(dim+1) + (j) )) 

double	d_eukl( long i1, long i2, long dim, double* x = NULL ) {

	double	result = 0.0, tmp;

	for( long i=0; i<dim; i++ )
	{
		if( i1 == -1 ) tmp = x[i];
		else tmp = opt_( i1, i );
		result += (tmp-opt_(i2,i)) * (tmp-opt_(i2,i));
	};
	return sqrt(result);
};

#undef opt_

//----------------------------------------------------------------------

void	calc_a( long optima, long dim, double min_r, double max_r ) {
// calculate coefficients in array a

	long	i, j;
	double	min_d;
	
	for( i=0; i < optima; ++i )
	{
		if( optima == 1 )
			min_d = max_r;
		else
			min_d = d_eukl( (i==0 ? 1 : 0), i, dim );
		for( j=0; j < optima; ++j )
			if( j != i && d_eukl( i, j, dim ) < min_d )
				min_d = d_eukl( i, j, dim ); 
		
		if( min_d < min_r*2 ) 
		{
			message( "min_r too big (testeval)",MSG_INFORMATION); 
			min_d == min_r*2; // bad ini data
		};
		min_d -= min_r;
		if( min_d > max_r ) min_d = max_r;
		a[i] = 2.0 / ( min_d * min_d );
	};
};

 
//-----------------------------------------------------------------------


double	field_st( long opt, long dim, double* x ) {
// 'field strength' with respect to optimum opt
// returns a value from [0,1]

	double	d = d_eukl( -1, opt, dim, x );
	
	if( d < 1.0 / sqrt( 2 * a[opt] ) )
		return 1.0 - a[opt] * d*d;
	else if( d < sqrt( 2.0 / a[opt] ) ) 
		return	a[opt]*(d-sqrt( 2/a[opt] )) * (d-sqrt( 2/a[opt] ));
	else
		return -d*d;
};



//-------------------------------------------------------------------------


double	NDimRealElement::Value() {
// INI.Fnum() can be used to select a funcion
// INI.Lowb(i) and INI.Upb(i) are lower and upper bounds on x[i]
// INI.Param(i) is the ith constant parameter, INI.ParamNum() is the
// number of available parameters
// dim is the dimension of space

// INI.Param(0) : min value of attraction area radius
// INI.Param(1) : max value of attraction area radius
// INI.Param(2) : number of local optima (max MaxLocalOptima (defined above))
// INI.Param(3 + i*(dim+1) + j ) : x_j in i-th local optima (from [0:1])
// INI.Param(3 + i*(dim+1) + dim ) : value of i-th local optima (positive)


	const long	optima = (long)INI.Param(2);
	double		d, result=-1e30, f_i;
	long count;
	long i,j,k,z,ii;
	
	
	// --- Variables para el retardo
	clock_t 		rel1,rel2;
	double 	res;
	time_t t1,t2;
	rel1=clock();
	
		
		
	//for(ii=0;ii<100000;ii++) sqrt(ii);
	
	rel2=clock();
	res=(double)(rel2-rel1)/(double)CLOCKS_PER_SEC;		
	
	// --- Variables para el retardo
	/*clock_t 		rel1,rel2;
	double 	res;
	time_t t1,t2;
	rel1=clock();
	*/
		
	//printf("Estoy aquí\n");	
	if( a[0] == 0.0 )
		// this runs only when eval first called
		calc_a( optima, dim, INI.Param(0), INI.Param(1) );
	
	for( i=0; i<optima; ++i )
	{
		f_i = INI.Param(3 + i*(dim+1) + dim );
		if( (d = field_st( i, dim, x ) * f_i) > result )
			result = d;
	};
	
	/*rel2=clock();
 	res=(double)(rel2-rel1)/(double)CLOCKS_PER_SEC;
	while (res< 0.00003) {
		rel2=clock();
		res=(double)(rel2-rel1)/(double)CLOCKS_PER_SEC;
 	};*/
	
	//The usleep() function delays program execution for the given number of micro_seconds. 
	//A microsecond is one millionth of a second.

	//usleep(10000); 

	
	return result;
};


