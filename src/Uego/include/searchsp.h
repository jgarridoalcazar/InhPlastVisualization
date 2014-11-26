#ifndef SEARCHSP_H
#define SEARCHSP_H

////////////////////////////////////////////////////////////
// $Id: searchsp.h,v 2.5 1998/03/17 23:14:52 jelasity Exp $
// searchsp.h
// declaration for abstract search space
////////////////////////////////////////////////////////////
// modification history:
////////////////////////////////////////////////////////////


//-------------------------------------------------------------------------
//---- abstract base class for different search spaces
//-------------------------------------------------------------------------


class SearchSpElement {

protected:

	// --- objective function value -----------------------------------
	double		value;
	virtual double	Value() = 0; // re-calculates obj. funct. value

	// --- error handling flag ----------------------------------------
	char	FailFlag;

public:

	virtual ~SearchSpElement() {};

	double		CurrValue() { return value; };
	virtual double	UpdateValue() { return ( value = Value() ); };

	// --- definition of distance in space -----------------------------
	virtual double	Distance( SearchSpElement*, SearchSpElement* s=NULL)=0;

	// ----- virtual functions for constructions (they work like new) --
	virtual SearchSpElement* RandNew() = 0;
		// random element from the whole space
	virtual SearchSpElement* RandNew( short, SearchSpElement* s=NULL ) = 0;
		// random element in area given by center and radius index
	virtual SearchSpElement* MutateNew( short, SearchSpElement* s=NULL )=0;
		// mutation on given element (may depend on radius index)
	virtual SearchSpElement* BetweenNew( SearchSpElement*,
					     SearchSpElement* s=NULL ) = 0;
		// element on the section connecting the given elements

	// ----- Improve *this using given spec. rad. index and maxevals ---
	virtual long	Optimize( short, long ); // returns funct. evals

	// ----- it substitues the "virtual assignement operator" -----------
	virtual void	UpdateFrom( SearchSpElement* ) = 0;

	// ----- virtual "static" functions for Configure -------------------
	virtual double	Diameter( Ini* ) = 0; // diam using bounds in ini
	virtual double	v( double ) = 0; // speed with given radius

	// ----- saves the element to stream --------------------------------
	virtual void	Save( FILE* ) = 0;
	virtual void	Save2( FILE* ) = 0;

	char	Fail() { return FailFlag; };
	virtual void	GetX( double* y)=0;

};

#endif

