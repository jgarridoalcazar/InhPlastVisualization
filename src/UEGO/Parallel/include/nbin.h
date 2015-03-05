#ifndef NBIN_H
#define NBIN_H

////////////////////////////////////////////////////////////
// $Id: nbin.h,v 2.5 1998/03/17 23:14:51 jelasity Exp $
// nbin.h
// declaration for n dimensional binary search space
////////////////////////////////////////////////////////////
// modification history:
// 	Jelasity 98 01 10: Mut and FlipBit()
////////////////////////////////////////////////////////////


//--------------------------------------------------------------------------
//---- base class for binary search spaces
//--------------------------------------------------------------------------


class BinaryElement : public SearchSpElement {

protected:

	virtual double	Value();

	static	long	dim;
	static	long	Mut; // working var. for MutateNew()
	unsigned long	*x; // the compressed space element

	void	FlipBit( long i ) { x[i/32] ^= (1UL << (i%32)); };


public:

	BinaryElement( long ); // creates an uninitialized 'empty' element
	virtual	~BinaryElement() { if( x != NULL ) delete x; };

	char	X( long i ) { // i: [0,dim-1], result is 0 or 1
		if( (x[i/32] & (1UL << (i%32))) != 0 ) return 1;
		else return 0; };

	virtual double	Distance( SearchSpElement*, SearchSpElement* s=NULL );

	// ----- implementations of SearchSpElement virtual functions
	virtual SearchSpElement* RandNew();
	virtual SearchSpElement* RandNew( short, SearchSpElement* s=NULL );
	virtual SearchSpElement* RandNewParal(){};
	virtual SearchSpElement* RandNewParal( short, SearchSpElement* s=NULL ){};
	

	virtual SearchSpElement* MutateNew( short, SearchSpElement* s=NULL );
	virtual SearchSpElement* BetweenNew( SearchSpElement*,
					     SearchSpElement* s=NULL );

	virtual void	UpdateFrom( SearchSpElement* );

	virtual double	Diameter( Ini* ini ) { return dim; };
	virtual double	v( double r ) { return 3.0/11.0 * sqrt(r); };
	virtual void	Save( FILE* );
	virtual void	Save2( FILE* ){};
	virtual void 	GetX( double* newx) {};

	virtual void SetX( double* newx) {};
	virtual void SetValue( double newx) {};
	virtual SearchSpElement* SetNew(double val,double* x){}; 
	virtual SearchSpElement* updateCenter(double *, double  ){};
	//virtual long	Optimize2( short, long ){}; // returns funct. evals


};


#endif

