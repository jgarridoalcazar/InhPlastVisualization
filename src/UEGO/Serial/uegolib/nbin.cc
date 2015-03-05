#include "uego.h"



long	BinaryElement::dim = 1;
long	BinaryElement::Mut = 0;

// -------------------------------------------------------------------------


BinaryElement::BinaryElement( long dimension ) {
// 'long' is assumed to be at least 32 bit long

	dim = dimension;

	FailFlag = ((x = new unsigned long[ (dim-1)/32 + 1 ]) == NULL);
	if( FailFlag ) message((char*)"No memory for BinaryElement.",MSG_ERROR);
};


// -------------------------------------------------------------------------


double	BinaryElement::Distance( SearchSpElement* e1, SearchSpElement* e2 ) {
// Hamming distance

        long		dist;
	BinaryElement	*r1, *r2;

	r1 = (BinaryElement*)e1;
	r2 = (BinaryElement*)e2;
	if( r2 == NULL ) r2 = this; // called with 1 par.

	dist = 0;
	for( long i=0; i < dim; ++i )
		if( r1->X(i) != r2->X(i) ) ++dist;

	return dist;
};


// -------------------------------------------------------------------------


SearchSpElement* BinaryElement::RandNew() {
// random bits

	long		i, j, k;
	BinaryElement	*result;

	FailFlag = 1==1;

	// --- new empty element
	result = new BinaryElement( dim );
	if( result == NULL || result->Fail() )
	{
		message( (char*)"Error creating empty element in BinaryElement::RandNew",MSG_ERROR);
		return (SearchSpElement*)NULL;
	};

	// --- filling the empty element
	j = 0;
	for( i=0; i <= (dim-1)/32; ++i )
	{
		result->x[i] = 0;
		for( k=0; k<32 && j<dim; ++k,++j )
		{
			result->x[i] <<= 1;
			result->x[i] |= (unsigned long) UegoLongRand( 0, 1 );
		};
	};
	result->UpdateValue();
	result->FailFlag = 1==0;

	FailFlag = 1==0;
	return result;
};


// -------------------------------------------------------------------------


SearchSpElement* BinaryElement::RandNew( short radind, SearchSpElement* e ) {
// max R(radind) number of bit mutations

	BinaryElement	*result, *r;
	long		i, bits;

	if( radind == 0 ) return RandNew(); // default for root

	FailFlag = 1==1;

	r = (BinaryElement*)e;
	if( r == NULL ) r = this; // called without par.

	// --- new empty element
	result = new BinaryElement( dim );
	if( result == NULL || result->Fail() )
	{
		message((char*) "Error creating empty element in BinaryElement::RandNew",MSG_ERROR);
		return (SearchSpElement*)NULL;
	};

	// --- filling the empty element
	result->UpdateFrom( r );
	bits = UegoLongRand( 1, ceil(INI.R( radind )) );
	for( i=0; i < bits; ++i ) // bit mutations
	{
		result->FlipBit( UegoLongRand( 0, dim-1 ) );
	};
	result->UpdateValue();
	result->FailFlag = 1==0;

	FailFlag = 1==0;
	return result;
};


// -------------------------------------------------------------------------


SearchSpElement* BinaryElement::MutateNew( short radind, SearchSpElement* e ) {
// mutation used by several algorithms (SHC, GAS, etc.)
// currently mutation probability (determined by MutJump) is fixed

	BinaryElement	*result, *r;
	const long	MutJump = dim/2; // mutation probability is fixed
	long		i;
	char		waschange = 1==0;

	FailFlag = 1==1;

	r = (BinaryElement*)e;
	if( r == NULL ) r = this; // called without par.

	// --- new empty element
	result = new BinaryElement( dim );
	if( result == NULL || result->Fail() )
	{
		message((char*) "Error creating empty element in BinaryElement::MutateNew",MSG_ERROR);
		return (SearchSpElement*)NULL;
	};

	// --- filling the empty element
	result->UpdateFrom( r );
	while( Mut < dim )
	{
		waschange = 1==1;
		result->FlipBit( Mut );
		Mut += UegoLongRand( 1, MutJump );
	};
	Mut -= dim;
	if( !waschange ) result->FlipBit( UegoLongRand( 0, dim-1 ) );

	result->UpdateValue();
	result->FailFlag = 1==0;

	FailFlag = 1==0;
	return result;
};


// -------------------------------------------------------------------------


SearchSpElement* BinaryElement::BetweenNew( SearchSpElement* e1, SearchSpElement* e2 ) {
// uniform crossover

	BinaryElement	*result, *r1, *r2;
	unsigned long	bit, i;

	FailFlag = 1==1;

	r1 = (BinaryElement*)e1;
	r2 = (BinaryElement*)e2;
	if( r2 == NULL ) r2 = this; // called with 1 par.

	// --- new empty element
	result = new BinaryElement( dim );
	if( result == NULL || result->Fail() )
	{
		message((char*) "Error creating empty element in BinaryElement::BetweenNew",MSG_ERROR);
		return (SearchSpElement*)NULL;
	};

	// --- filling the empty element
	for( i=0; i <= (dim-1)/32; ++i ) result->x[i] = 0;
	for( i=0; i < dim; ++i )
	{
		if( UegoLongRand(0,1) == 0 ) bit = r1->X(i);
		else bit = r2->X(i);
		result->x[i/32] |= bit << (i%32);
	};
	result->UpdateValue();
	result->FailFlag = 1==0;

	FailFlag = 1==0;
	return result;
};

// -------------------------------------------------------------------------


void	BinaryElement::UpdateFrom( SearchSpElement* e ) {

	value = e->CurrValue();
	for( long i=0; i <= (dim-1)/32; ++i )
	{
		x[i] = ( (BinaryElement*)e )->x[i];
	};
};


// -------------------------------------------------------------------------


void	BinaryElement::Save( FILE* stream ) {

	FailFlag = 1==0;

	for( long i=dim-1; i>=0; --i )
	{
		FailFlag |= fprintf( stream, "%ld", (long)X(i) ) == EOF;
	};

	if( value < 1e10 )
		FailFlag |= fprintf( stream, " %.10le\n", value ) == EOF;
	else
		FailFlag |= fprintf( stream, " %lf\n", value ) == EOF;
	if( Fail() ) message((char*)"Error writing BinaryElement.",MSG_ERROR);
};

