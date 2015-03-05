#include "uego.h"

double	BinaryElement::Value() {
// INI.Fnum() can be used to select a funcion
// INI.Param(i) is the ith constant parameter, INI.ParamNum() is the
// number of available parameters
// X(i) is the i-th bit, 0 or 1; dim is the length of bitstring

	long	i, y;
	double  n1, n2;

	switch( INI.Fnum() ) {

	case 0: // subset sum
		n1 = INI.Param(0); // sum to create
		n2 = 0;
		for( i=0; i<dim; ++i )
			if( X(i)==1 ) n2 += INI.Param(i+1);
		if( n2 <= n1 ) return n2-n1;
		else return -n2;

	default: // bit counting
		y = 0;
		for( i=0; i<dim; ++i )
			if( X(i) == 1 ) ++y;
		return y;
	};
};
