#include "uego.h"

#include <fstream>

long	NDimRealElement::dim = 1;

// -------------------------------------------------------------------------


NDimRealElement::NDimRealElement( long dimension ) {

	b = NULL;
	epsi = NULL;

	dim = dimension;

	FailFlag = ((x = new double[ dim * 2 ]) == NULL);
	if( FailFlag ) message((char*)"No memory for NDimRealElement.",MSG_ERROR);
};


// -------------------------------------------------------------------------


double	NDimRealElement::UpdateValue() {

	double	*normx = x;

	// --- converting to real coordinates
	x += dim;
	for( long i=0; i < dim; ++i )
	{
		x[i] = normx[i]*(INI.Upb(i)-INI.Lowb(i)) + INI.Lowb(i);
	};

	// --- counting value
	value = Value();

	// --- restoring normalized coordinates
	x = normx;

	return value;
};


// -------------------------------------------------------------------------


double	NDimRealElement::Diameter( Ini* ini ) {
// Euclidian distance of lower and upper bounds

	double	sqrsum = 0.0;

	for( long i=0; i < dim; ++i )
	{
		sqrsum += ( ini->Lowb(i) - ini->Upb(i) ) *
				( ini->Lowb(i) - ini->Upb(i) );
	};

	return sqrt( sqrsum ); 
};


// -------------------------------------------------------------------------


double	NDimRealElement::v( double r ) {
#define SPEED_CONSTANTS 14
//  Using that binom(x,(x-1)/2)==.5*binom(x+1,(x+1)/2) and
//  binom(2n,n) is approx. 2**(2n)/sqrt(pi*n)

	double	V[ SPEED_CONSTANTS ] = {	// for r==1.0
		0.0,			.25,
		M_2_PI / 3.0,		3.0 / 16.0,
		4 * M_2_PI / 15.0,	10.0 / 64.0,
		8 * M_2_PI / 35.0,	35.0 / 256.0,
		64 * M_2_PI / 315.0,	63.0 / 512.0,
		128 * M_2_PI / 693.0,	231.0 / 2048.0,
		512 * M_2_PI / 3003.0,	858.0 / 8192.0  };

	if( dim < SPEED_CONSTANTS ) return V[dim] * r;
	else return r / sqrt( M_PI * (2 * (dim+1)) );

#undef SPEED_CONSTANTS
};


// -------------------------------------------------------------------------


double	NDimRealElement::Distance( SearchSpElement* e1, SearchSpElement* e2 ) {
// Euclidian distance in real coordinate space

	double		sqrsum, d1, d2;
	NDimRealElement *r1, *r2;

	r1 = (NDimRealElement*)e1;
	r2 = (NDimRealElement*)e2;
	if( r2 == NULL ) r2 = this; // called with 1 par.

	sqrsum = 0.0;
	for( long i=0; i < dim; ++i )
	{
		//d1 = r1->x[i]*(INI.Upb(i)-INI.Lowb(i)) + INI.Lowb(i);
		d1 = r1->x[i];
		//d2 = r2->x[i]*(INI.Upb(i)-INI.Lowb(i)) + INI.Lowb(i);
		d2 = r2->x[i];
		sqrsum += ( d1 - d2 ) * ( d1 - d2 );
	};

	return sqrt( sqrsum );
};


// -------------------------------------------------------------------------


SearchSpElement* NDimRealElement::RandNew() {
// uniform distribution in search space

	NDimRealElement	*result;
	long flipBit;

	FailFlag = 1==1;

	// --- new empty element
	result = new NDimRealElement( dim );
	if( result == NULL || result->Fail() )
	{
		message( (char*)"No memory in NDimRealElement::RandNew(1)",MSG_ERROR);
		return (SearchSpElement*)NULL;
	};

	// --- Initially, we consider all the spikes  ares considered.
	// --- Filling the empty element

	for( long i=0; i < dim; ++i )
		result->x[i] = UegoRand() ; //drand48(); 
	
	result->UpdateValue();
	
	FailFlag = 1==0;
	return result;
};


// -------------------------------------------------------------------------

SearchSpElement* NDimRealElement::RandNewParal() {
// uniform distribution in search space

	NDimRealElement	*result;
	long flipBit;

	FailFlag = 1==1;

	// --- new empty element
	result = new NDimRealElement( dim );
	if( result == NULL || result->Fail() )
	{
		message( (char*)"No memory in NDimRealElement::RandNew(1)",MSG_ERROR);
		return (SearchSpElement*)NULL;
	};

	// --- Initially, we consider all the spikes  ares considered.
	// --- Filling the empty element

	for( long i=0; i < dim; ++i )
		result->x[i] = UegoRand() ; //drand48(); 
	
	//result->UpdateValue();
	
	FailFlag = 1==0;
	return result;
};


// -------------------------------------------------------------------------


SearchSpElement* NDimRealElement::RandNew( short radind, SearchSpElement* e ) {
// uniform distribution in given area (a hyper-rectangle)
	NDimRealElement *result, *r;
	double		rad;
	long		flipBit;

	if( radind == 0 ) return RandNew(); // default for root
	
	rad = INI.R( radind ) / INI.R( 0 ) * sqrt( dim ); // normalized radius

	FailFlag = 1==1;

	r = (NDimRealElement*)e;
	if( r == NULL ) r = this; // called without par.

	// --- new empty element
	result = new NDimRealElement( dim );
	if( result == NULL || result->Fail() )
	{
		message((char*) "No memory in NDimRealElement::RandNew(2)",MSG_ERROR);
		return (SearchSpElement*)NULL;
	};

	// --- Filling the empty element
	for( long i = 0; i < dim; ++i )
	{
		result->x[i] = UegoDoubleRand( r->x[i]-rad, r->x[i]+rad );
		if( result->x[i] < 0.0 ) result->x[i] = 0.0;
		else if( result->x[i] > 1.0 ) result->x[i] = 1.0;
	};

	result->UpdateValue();

	FailFlag = 1==0;
	return result;
};

// -------------------------------------------------------------------------


SearchSpElement* NDimRealElement::RandNewParal( short radind, SearchSpElement* e ) {
// uniform distribution in given area (a hyper-rectangle)
	NDimRealElement *result, *r;
	double		rad;
	long		flipBit;

	if( radind == 0 ) return RandNewParal(); // default for root
	
	//rad = INI.R( radind ) / INI.R( 0 ) * sqrt( dim ); // normalized radius
	rad = INI.R( radind );

	printf("Generating new species with rad=%e\n",rad);

	FailFlag = 1==1;

	r = (NDimRealElement*)e;
	if( r == NULL ) r = this; // called without par.

	// --- new empty element
	result = new NDimRealElement( dim );
	if( result == NULL || result->Fail() )
	{
		message((char*) "No memory in NDimRealElement::RandNew(2)",MSG_ERROR);
		return (SearchSpElement*)NULL;
	};

	// --- Filling the empty element
	for( long i = 0; i < dim; ++i )
	{
		result->x[i] = UegoDoubleRand( r->x[i]-rad, r->x[i]+rad );
		if( result->x[i] < 0.0 ) result->x[i] = 0.0;
		else if( result->x[i] > 1.0 ) result->x[i] = 1.0;
	};

	//result->UpdateValue();

	FailFlag = 1==0;
	return result;
};


// -------------------------------------------------------------------------

/*
void	NDimRealElement::SetNew( double* y) {
			
	for( long i=0; i < dim; ++i )
	{
		x[i] += sign * y[i];
		if( x[i] < 0.0 ) x[i] = 0.0;
		else if( x[i] > 1.0 ) x[i] = 1.0;
	};
}*/

// -------------------------------------------------------------------------


/*SearchSpElement* NDimRealElement::SetNew(double val ,double* xx) {
// set the values of the new SearchSpElement

	NDimRealElement *result, *r;
	double		rad;

	
	FailFlag = 1==1;


	// --- new empty element
	result = new NDimRealElement( dim );
	if( result == NULL || result->Fail() )
	{
		message( "No memory in NDimRealElement::SetNew(2)",MSG_ERROR);
		return (SearchSpElement*)NULL;
	};

	// --- filling the empty element
	for( long i=0; i < dim; ++i )
	{
		result->x[i] = xx[i];

	
	};
	result->SetValue(val);
	//El level se inicializa en la funcion NewSpecies del master. Dentro de : newevals = newspechead.NewSpecies( &newlst, tmplevel,level, tempeval, &sizelist );
					
	FailFlag = 1==0;
	return result;
};
*/


SearchSpElement* NDimRealElement::SetNew( double fobjValue, double* xx) {
// set the values of the new SearchSpElement

	NDimRealElement *result=NULL;
	double		rad;
	
	FailFlag = 1==1;

	// --- new empty element
	result = new NDimRealElement( dim );	
	if( result == NULL || result->Fail() )
	{
		message((char*) "No memory in NDimRealElement::SetNew(2)",MSG_ERROR);
		return (NDimRealElement*)NULL;
	};

	// --- filling the empty element
	for( long i=0; i < dim; ++i )
		result->x[i] = xx[i];  
	
	//-- actualizamos el valor de la funcion objetivo
	result->value = fobjValue;	
	this->UpdateFrom(result);
	
	//printf("SetNew == (%f, %f, %f) = %f\n", x[0], x[1], x[2], value);
	
	FailFlag = 1==0;
	return result;
};

// -------------------------------------------------------------------------

SearchSpElement* NDimRealElement::updateCenter(double *aux, double fobjValue ) {
// uniform distribution in search space

	NDimRealElement	*result=NULL;

	FailFlag = 1==1;

	// --- new empty element
	result = new NDimRealElement( dim );
	if( result == NULL || result->Fail() )
	{
		message((char*) "No memory in NDimRealElement::updateCenter()",MSG_ERROR);
		return (SearchSpElement*)NULL;
	};

	// --- filling the empty element	
	for( long i=0; i < dim; ++i )
		result->x[i] =aux[i]; 
	
	//-- actualizamos el valor de la funcion objetivo
	result->value = fobjValue;	
	//printf("||||||||||value %lf fobjValue %lf\n", value, fobjValue);

	FailFlag = 1==0;
	return result;
};

// -------------------------------------------------------------------------


void	NDimRealElement::Add( double* y, signed char sign ) {
			
	for( long i=0; i < dim; ++i )
	{
		x[i] += sign * y[i];
		if( x[i] < 0.0 ) x[i] = 0.0;
		else if( x[i] > 1.0 ) x[i] = 1.0;
	};
}


// -------------------------------------------------------------------------


double 	NDimRealElement::Gauss( double bias, double sigma ) { 

	double	xx,
		z = 0.0;

	for( int k=0; k<12; k++ ) 
	{
	   xx = UegoRand(); 
      	   z = z + xx;
	}

	return ( bias + sigma * (z-6.0) );
};



// -------------------------------------------------------------------------


SearchSpElement* NDimRealElement::MutateNew( short radind, SearchSpElement* e ) {
// mutation used by several algorithms (SHC, GAS, etc.)

return RandNew( radind, e );
// ok, just to make it work quickly, temporarily
};


// -------------------------------------------------------------------------


SearchSpElement* NDimRealElement::MutateNewParal( short radind, SearchSpElement* e ) {
// mutation used by several algorithms (SHC, GAS, etc.)

return RandNewParal( radind, e );
// ok, just to make it work quickly, temporarily
};


// -------------------------------------------------------------------------


SearchSpElement* NDimRealElement::BetweenNew( SearchSpElement* e1, SearchSpElement* e2 ) {

	NDimRealElement *result, *r1, *r2;
	
	FailFlag = 1==1;

	r1 = (NDimRealElement*)e1;
	r2 = (NDimRealElement*)e2;
	if( r2 == NULL ) r2 = this; // called with 1 par.

	// --- new empty element
	result = new NDimRealElement( dim );
	if( result == NULL || result->Fail() )
	{
		message((char*) "No memory in NDimRealElement::BetweenNew",MSG_ERROR);
		return (SearchSpElement*)NULL;
	};

	
	// --- filling the empty element

	for( long i=0; i < dim; ++i )
		result->x[i] = ( r1->x[i] + r2->x[i] ) / 2.0;
	result->UpdateValue();

	FailFlag = 1==0;
	return result;
};


// -------------------------------------------------------------------------


void	NDimRealElement::UpdateFrom( SearchSpElement* e ) {

	value = e->CurrValue();
	for( long i=0; i < dim; ++i )
	{
		x[i] = ( (NDimRealElement*)e )->x[i];
	};
};


// -------------------------------------------------------------------------
void	NDimRealElement::Save( FILE* stream ) {

	double	y;

	FailFlag = 1==0;

	// --- saving as real coordinates
	for( long i=0; i < INI.Dimension(); ++i )
	{
		y =x[i]*(INI.Upb(i)-INI.Lowb(i)) + INI.Lowb(i);
		FailFlag |= fprintf( stream, "%.10le ", y ) == EOF;
	};

	FailFlag |= fprintf( stream, "%.10le\n", value ) == EOF;
	if( Fail() ) message((char*)"Error writing SearchSpElement.",MSG_ERROR);
};


// -------------------------------------------------------------------------
std::ofstream&	NDimRealElement::Save( std::ofstream & myfile ) {

	double	y;

	// --- saving as real coordinates
	for (unsigned int i=0; i<INI.Dimension(); ++i){
		if (INI.ParameterScale(i)==LOGARITHMIC){
			double logmin = log10(INI.Lowb(i));
			double logmax = log10(INI.Upb(i));
			y= pow(10.0,(x[i]*(logmax - logmin)))*INI.Lowb(i);
		} else {
			y=x[i]*(INI.Upb(i)-INI.Lowb(i)) + INI.Lowb(i);
		}
		myfile << y << "\t";
	}

	myfile << value;

	if (!myfile){
		message((char*)"Error writing SearchSpElement.",MSG_ERROR);
	}

	return myfile;

};

// -------------------------------------------------------------------------
SearchSpElement *	NDimRealElement::LoadFromFile(std::ifstream & file){

	double y;

	NDimRealElement * NewElement = new NDimRealElement(INI.Dimension());

	for (long i=0; i < INI.Dimension(); ++i ){
		file >> y;

		if (!file){
			message((char*)"Error reading SearchSpElement.",MSG_ERROR);
		}

		// Normalize the values
		if (INI.ParameterScale(i)==LOGARITHMIC){
			NewElement->x[i] = log10(y/INI.Lowb(i))/log10(INI.Upb(i)/INI.Lowb(i));
		} else {
			NewElement->x[i] = (y-INI.Lowb(i))/(INI.Upb(i)-INI.Lowb(i));
		}

	}

	file >> NewElement->value;

	if (!file){
		message((char*)"Error reading SearchSpElement.",MSG_ERROR);
	}

	return NewElement;
}

// -------------------------------------------------------------------------
void	NDimRealElement::Save2( FILE* stream ) {
	
	double	y;
	int 	yaux;
	FailFlag = 1==0;

	long numActiveRBFs;
	double minSigma, avSigma,maxSigma,DevSigma;
	long 	aux;
		

	numActiveRBFs=0;

	aux = dim/3;
	
	// --- saving as real coordinates
	//---Number of active RBFs
	for( long i=0; i < aux; ++i )
	{	
		yaux=(int)(x[i]+0.5);
		if(yaux==1) numActiveRBFs ++;
		FailFlag |= fprintf( stream, "%i  ", yaux ) == EOF;
	}
	//---Sigmas
	//fprintf(stream,"\n");
	minSigma=INFINITY;
	maxSigma=-INFINITY;
	avSigma=0;
	DevSigma=0;

	for( long i=aux; i < (2*aux); ++i )
	{
		y = x[i]*(INI.Upb(i)-INI.Lowb(i)) + INI.Lowb(i);
		if(y < minSigma) minSigma = y;
		if(y > maxSigma) maxSigma = y;
		avSigma += y;		
 
		FailFlag |= fprintf( stream, "%.10le  ", y ) == EOF;
	};

	avSigma/= aux; 


	//---Centers
	for( long i=(2*aux); i < dim; ++i )
	{
		y = x[i]*(INI.Upb(i)-INI.Lowb(i)) + INI.Lowb(i);
		FailFlag |= fprintf( stream, "%.10le  ", y ) == EOF;
	};

	//fprintf(stream,"\n");

	//---Compute DevSigma
	for( long i=aux; i < (2*aux); ++i )
	{
		y = x[i]*(INI.Upb(i)-INI.Lowb(i)) + INI.Lowb(i);
		DevSigma+= ( (y-avSigma) * (y-avSigma));		
 
	};
	DevSigma/=aux;
	DevSigma = sqrt(DevSigma);

	fprintf(stream,"\n");
	//---Objective function
	FailFlag |= fprintf( stream, "%.10le\n", value ) == EOF;
	if( Fail() ) message((char*) "Error writing NDimRealElement.",MSG_ERROR);


	//---Metrics
	FailFlag |= fprintf( stream, "%ld %.10le %.10le %.10le %.10le\n \n", numActiveRBFs,minSigma,avSigma,maxSigma,DevSigma ) == EOF;
	if( Fail() ) message((char*) "Error writing NDimRealElement.",MSG_ERROR);



};

