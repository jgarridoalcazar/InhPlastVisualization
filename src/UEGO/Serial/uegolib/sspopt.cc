#include "uego.h"


long	SearchSpElement::Optimize( short radind, long maxeval ) {

	SearchSpElement	*spe;

	FailFlag = 1==1;

	for( long i=0; i<maxeval; ++i )
	{
		spe = MutateNew( radind );
		if( Fail() ) return -1;
		else FailFlag = 1==1;

		if( spe->CurrValue() >= CurrValue() ) UpdateFrom( spe );
		delete spe;
	};

	FailFlag = 1==0;

	return maxeval;
};

