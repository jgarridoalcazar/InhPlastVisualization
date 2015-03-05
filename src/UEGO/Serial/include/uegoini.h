#ifndef INI_H
#define INI_H

////////////////////////////////////////////////////////////
// $Id: uegoini.h,v 2.5 1998/03/17 23:14:52 jelasity Exp $
// uegoini.h
// contains the declaration of class Ini
// Ini contains and handles the settings used in an uego
// session.
////////////////////////////////////////////////////////////
// modification history:
//	jelasity 98 02 20 name changed to uegoini.h
////////////////////////////////////////////////////////////


class Ini {
friend class Configure;

private:

	// --- random numbers --------------------------------------
	unsigned long	seed;	// random seed

	// --- objective function -------------------------------------
	long	type;		// type of search space
	long	fnum;		// number of obj. funtc. in library
	long	dimension;	// dimesion of the search space
	long	paramnum;	// length of param[]
	double	*param;		// optional parameters fot the obj. funct.
	double	*lowb;		// lower bounds of variables for real spaces
	double	*upb;		// upper bounds of variables for real spaces
	char    **paramNames;	// parameter names
	SearchSpElement	*prototype;	// points to an instance of the
					// class selected by the user

	// --- uego values -----------------------------------------------
	unsigned long	maxevals;	// max num. of funct. evals
	long	maxspecnum;		// max number of species at a time
	double	threshold;		// stability threshold
	short	levels;			// max strict level (length of r[])

	// --- automatic
	double	*r;			// radii of different levels
	unsigned long	*evals;		// max funct. evals on  different levels
	unsigned long	*newspecevals;	// max funct. evals in NewSpecies

	// --- handy things for Save and Ini( FILE* ) ----------------------
	static char	SaveVector( FILE*, char*, double*, unsigned long*, long );
	static char	ToBuffer( FILE*, char** );
	static char	GetVector( char*, char*, double**, unsigned long**, long );
	double	GetValue( char*, char* );
	char**    GetNames(char*, char*,long);
	char	FailFlag;

public:

	Ini( FILE* );
	Ini();
	~Ini();

	unsigned long	*Seed() { return &seed; };
		// random seed; modified by UegoRand()
		// that's why it's a pointer

	long	Type() { return type; };
	long	Fnum() { return fnum; };
	long	Dimension() { return dimension; };
	long	ParamNum() { return paramnum; };
	double	Param( long i ) {
		if( i < 0 ) i = 0; else if( i>=paramnum ) i=paramnum-1;
		return param == NULL ? 0.0 : param[i]; };
	double	Lowb( long i ) {
		if( i < 0 ) i = 0; else if( i>=dimension ) i=dimension-1;
		return lowb == NULL ? 0.0 : lowb[i]; };
	double	Upb( long i ) {
		if( i < 0 ) i = 0; else if( i>=dimension ) i=dimension-1;
		return upb == NULL ? 0.0 : upb[i]; };

	char * ParameterName(long i) {
		if( i < 0 ) i = 0; else if( i>=dimension ) i=dimension-1;
		return paramNames[i]; };


	unsigned long	MaxEvals() { return maxevals; };
	long		MaxSpecNumber() { return maxspecnum; };
	short		Levels() { return levels; };
	double		R( long i ) {
		if( i < 0 ) i = 0; else if( i>=levels ) i=levels-1;
		return r == NULL ? 0.0 : r[i]; };
	unsigned long	Evals( long i ) {
		if( i < 0 ) i = 0; else if( i>=levels ) i=levels-1;
		return evals == NULL ? 0 : evals[i]; };
	unsigned long	NewSpecEvals( long i ) {
		if( i < 1 ) i = 1; else if( i>=levels ) i=levels-1;
		return newspecevals == NULL ? 0 : newspecevals[i-1]; };

	SearchSpElement	*Prototype() { return prototype; };
	void		SetPrototype(); // creates prototype of type 'type'

	void	Save( FILE* );

	char	Fail() { return FailFlag; };
	
	double	v( double r );
};

#endif
