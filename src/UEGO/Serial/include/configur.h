#ifndef CONFIGUR_H
#define CONFIGUR_H

////////////////////////////////////////////////////////////
// $Id: configur.h,v 2.5 1998/03/17 23:14:51 jelasity Exp $
// configur.h
// contains the declaration of class Configure
// Configure is completely static and is used for creating
// uego configuration files using an ini object as starting
// point and for saving the result. It is separated from
// class Ini for clarity only.
////////////////////////////////////////////////////////////
// modification history:
//	jelasity 98 01 17 command line interface
//	jelasity 98 01 23 ReAllocVectors has a long par.
////////////////////////////////////////////////////////////

#define NOTHING_MISSING	-1	// maxevals
#define MISSING_N	0	// maxevals
#define MISSING_M	1	// max number of species
#define MISSING_NU	2	// threshold
#define MISSING_L	3	// levels
#define MISSING_LAST_R	4	// smallest radius
#define VALID_MISS(x) ( (x)>=0 && (x)<=4 && (x)!=3 )

class Configure {

private:
	// --- first default then result ini
	static Ini	*ini;

	// --- I/O things -----------------------------------------
	static double	GetValue( char*, char*, double );
	static FILE*	GetOpenFile( char*, char* );
	static void	ReadParFile( FILE* );
	static void	GetParameterVector();
	static void	GetBoundVectors( long );

	// --- automatic parameter setting things
	static double	last_r;
	static long	Diff();
	static void	SetMissingN();
	static void	SetMissingNu();
	static void	SetMissingLastR();
	static void	SetMissingL();
	static void	SetMissingM();
	static void	SetVectors( long );
	static void	ReAllocVectors( long );
	static void	AutoParameters( long );

	// --- configuring things
	static void	DisplayConfig();
	static void	ParamConfig( int, char** );

	static char	FailFlag;

public:

	static void	Config( char*, int, char** );

	static char	Fail() { return FailFlag; };
};


#endif
