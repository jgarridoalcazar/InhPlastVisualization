#include "uego.h"

long	NDimRealElement::dim = 1;

// -------------------------------------------------------------------------


NDimRealElement::NDimRealElement( long dimension ) {

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
		d1 = r1->x[i]*(INI.Upb(i)-INI.Lowb(i)) + INI.Lowb(i);
		d2 = r2->x[i]*(INI.Upb(i)-INI.Lowb(i)) + INI.Lowb(i);
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
	for( long i = 0; i < dim/3; ++i ) 
		result->x[i] = 1.0;
	
	// --- Randomly select a spike to be zero
	if(UegoRand() > 0.5){
		flipBit = UegoLongRand(0, ((dim/3)-1));	
		result->x[flipBit] = 0.0;
	}

	// --- filling the empty element
	for( long i = dim/3; i < dim; ++i ) 
		result->x[i] = UegoRand();
	
	result->UpdateValue();
	
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
	
	//rad = INI.R( radind ) / INI.R( 0 ) * sqrt( dim ); // normalized radius

	//-------------------------------------------
	//----NORMALIZAMOS RADIO!!!!!!!!!!!!!!!!!!!!!
	//-------------------------------------------
	rad=INI.R(radind)/INI.R(0)*sqrt((dim/3)*2), //normalized rad.
			


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
	for( long i = 0; i < dim/3; ++i ) 
		result->x[i] = 1.0;
	
	// --- Randomly select a spike to be zero
	if(UegoRand() > 0.5){
		flipBit = UegoLongRand(0, ((dim/3)-1));	
		result->x[flipBit] = 0.0;
	}
	for( long i = dim/3; i < dim; ++i )
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


void	NDimRealElement::Add( double* y, signed char sign ) {
			
	for( long i=dim/3; i < dim; ++i )
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

	for( long i=0; i < dim/3; ++i ){
		if( result->x[i] == 0.5){
			if(UegoRand() > 0.5) result->x[i] =1;
			else result->x[i] =0;
		}			 		
	}
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

