#ifndef SPECLIST_H
#define SPECLIST_H

////////////////////////////////////////////////////////////
// $Id: speclist.h,v 2.5 1998/03/17 23:14:52 jelasity Exp $
// speclist.h
// contains the declaration of class SpeciesList
////////////////////////////////////////////////////////////
// modification history:
//
////////////////////////////////////////////////////////////

class SpeciesList {
friend class Master;

private:

	SearchSpElement	*center;	// center of species
	short		level;		// strict level of species

	// ----- methods for the main algorithm: Optimize and NewSpecies
	// ----- value of center may be changed in both
	long	Optimize(  long = -1 );	// -1 : use Evals()
	long	Optimize2(  long = -1 );	// -1 : use Evals()
	
	long	Evals( long = -1 );// evaluations  (depends on level and ini)
	long	_NewSpecies( SpeciesList**, short, long );
		// returns list of new species that have the given strict level
		// using maximum given function evaluations
	long	NewSpecies( SpeciesList**, short, long );

	SpeciesList	*prev;
	SpeciesList	*next;

	char	FailFlag;

public:

	SpeciesList() { prev = next = NULL; center = NULL; FailFlag = 1==0; };
	SpeciesList( SearchSpElement* c, short l ) : center(c), level(l) {
		prev = next = NULL; FailFlag = 1==0; };
	~SpeciesList() { if( center != NULL ) delete center; };

	void	Save( FILE* );
	void	Save2( FILE* );
	void	Absorb( SpeciesList* ); // absorb given species

	char	Fail() { return FailFlag; };
};


#endif
