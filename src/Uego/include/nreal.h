#ifndef NREAL_H
#define NREAL_H

////////////////////////////////////////////////////////////
// $Id: nreal.h,v 2.5 1998/03/17 23:14:52 jelasity Exp $
// nreal.h
// declaration for n dimensional real search spaces
////////////////////////////////////////////////////////////
// modification history:
//	Jelasity 98 02 15 inner repr. is (0,1)^dim
////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------
//---- base class for n dim real search spaces
//--------------------------------------------------------------------------
// the inner representation works on a (0,1)^n cube, but Value() sees
// the real user-given values in x. This is implemented through UpdateValue().
// Value() is never called directly. The inner representation does not affect
// the behaviour of public functions.
//--------------------------------------------------------------------------

class NDimRealElement : public SearchSpElement {

protected:

	virtual double	Value();
	double		Gauss( double, double );
	void		Add( double*, signed char = 1 );
	static long	dim;
	double		*x; // the space element

public:

	NDimRealElement( long ); // creates an uninitialized 'empty' element
	virtual ~NDimRealElement() { if( x != NULL ) delete x; };

	virtual double	UpdateValue();
	
	virtual double	Distance( SearchSpElement*, SearchSpElement* s=NULL );

	// ----- implementations of SearchSpElement virtual functions
	virtual SearchSpElement* RandNew();
	virtual SearchSpElement* RandNew( short, SearchSpElement* s=NULL );
	virtual SearchSpElement* MutateNew( short, SearchSpElement* s=NULL );
	virtual SearchSpElement* BetweenNew( SearchSpElement*,
					     SearchSpElement* s=NULL );

	virtual long	Optimize( short, long ); // returns funct. evals
	
	virtual void	UpdateFrom( SearchSpElement* );

	virtual double	Diameter( Ini* );
	virtual double	v( double );

	virtual void	Save( FILE* );
	virtual void	Save2( FILE* );
	virtual void	GetX( double* y){for (int k=0; k<dim;k++) y[k]=x[k];};

};

#endif

