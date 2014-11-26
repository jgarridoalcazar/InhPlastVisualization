#include <string.h>
#include "uegoconf.h"
#include "uego.h"
#include "configur.h"



Ini*	Configure::ini;
char	Configure::FailFlag;
double	Configure::last_r;
char	prompt[101];


// -----------------------------------------------------------------------
// --- I/O functions to help reading consol and files
// -----------------------------------------------------------------------


double	Configure::GetValue( char* prompt, char* helpmsg, double def ) {
// purpose:
// 	read a real value from the standard input stream interactively
//	'prompt' is displayed as prompt for input request
//	if user gives h<enter>, shows 'helpmsg' and 'prompt' again
//	if user gives <enter> returns 'def'
//	if user gives q<enter> returns 0.0 and sets FailFlag.
// return:
//	on success: the value read, and FailFlag is set to false
//	on error (i.e. cancellation): FailFlag is set to true

	char	puffer[102], ValueRead;
	double	result;

	if( Fail() ) return 0.0;
	else FailFlag = 1==1;

	ValueRead = 1==0;
	while( !ValueRead )
	{
		printf( prompt );
		fgets( puffer, 100, stdin );
		switch( puffer[0] )
		{
			case 'q': // cancel config
				message("Configuration cancelled.",
							MSG_INFORMATION);
				return 0.0; // FailFlag is set to true
			case 'h': // show help message
				printf( helpmsg );
				break;
			case '\n': // accept def
				result = def;
				ValueRead = 1==1;
				break;
			default: // try to read value
				ValueRead =
				sscanf( puffer, "%lf" , &result ) == 1;
		};
	};

	FailFlag = 1==0;
	return result;
};


// -----------------------------------------------------------------------


FILE*	Configure::GetOpenFile( char* prompt, char* helpmsg ) {
// purpose:
// 	open a file using a name read from the standard input stream
//	interactively. 'prompt' is displayed as prompt for input request
//	if user gives h<enter>, shows 'helpmsg' and 'prompt' again
//	if user gives <enter> returns NULL;
//	if user gives q<enter> returns NULL and sets FailFlag.
// return:
//	on success: the opened file, and FailFlag is set to false
//	on error (i.e. cancellation): FailFlag is set to true

	char	puffer[102], ValueRead;
	FILE	*result;

	if( Fail() ) return NULL;
	else FailFlag = 1==1;

	ValueRead = 1==0;
	while( !ValueRead )
	{
		printf( prompt );
		fgets( puffer, 100, stdin );
		switch( puffer[0] )
		{
		case 'h': // show help message
				printf( helpmsg );
				break;
		case '\n': // accept def
				result = NULL;
				ValueRead = 1==1;
				break;
		case 'q': // cancel config if only q given
			if( puffer[1] == '\n' )
			{
				message("Configuration cancelled.",
							MSG_INFORMATION);
				return NULL; // FailFlag is set to true
			};
		default: // try to read value
			puffer[ strlen(puffer)-1 ] = 0; // clear '\n'
			ValueRead =
				(result = fopen( puffer, "rt" )) != NULL;
		};
	};

	FailFlag = 1==0;
	return result;
};


// -----------------------------------------------------------------------


void	Configure::ReadParFile( FILE* parfile ) {

	char	puffer[10002], ValueRead, *pos,
		white[15] = " \n\v\b\r\f\t";
	long	i;
	double	d;

	if( Fail() ) return;
	else FailFlag = 1==1;

	// --- get length --------------------------------------
	// ---read length; after reading pos shows the next part
	// ---of the line of length
	fgets( puffer, 10000, parfile );
	while( !feof( parfile ) && puffer[0]=='#' )
		fgets( puffer, 10000, parfile );
	pos = puffer;
	if( !feof(parfile) )
		pos += strspn( pos, white ); // skip leading whitespace
	ValueRead = !feof(parfile) && sscanf( pos, "%ld" , &i ) == 1;
	if( ValueRead && i > 0 ) ini->paramnum = i;
	else
	{
		message( "Bad parameter file (length).",MSG_ERROR);
		return;
	};
	pos = strpbrk( pos, white ); // skip number just read


	// --- memory allocation -------------------------------
	ini->param = new double[ ini->paramnum ];
	if( ini->param == NULL )
	{
		message( "No memory in configuration.",MSG_ERROR);
		return;
	};

	// --- read parameters --------------------------------------
	// (line of paramnum in puffer here)
	i = 0;
	while( !feof(parfile) && i<ini->paramnum )
	{
		pos += strspn( pos, white ); // skip leading whitespace
		ValueRead = sscanf( pos, "%lf ", &d ) == 1;
		pos = strpbrk( pos, white ); // skip number just read
		if( !ValueRead )
		{
			// --- read next line
			fgets( puffer, 10000, parfile );
			while( !feof( parfile ) && puffer[0]=='#' )
				fgets( puffer, 10000, parfile );
			pos = puffer;
		}
		else ini->param[ i++ ] = d;
	};
	if( i < ini->paramnum )
	{
		message("Bad parameter file (parameters).",MSG_ERROR);
		return;
	};

	FailFlag = 1==0;
};


// -----------------------------------------------------------------------


void	Configure::GetParameterVector() {

	FILE	*parfile;
	long	newparnum, i;

	if( Fail() ) return;

	// --- trying to get a filename to open --------------------------
	parfile = GetOpenFile( "Parameter text file? (<no file>):", PARS_HE );
	if( Fail() ) return; // configure cancelled
	if( parfile != NULL )
	{
		ReadParFile( parfile );
		fclose( parfile );
		return;
	};

	// --- No filename, get number of pars ---------------------------
	sprintf( prompt, "%s (%ld): ", PARNUM_T, ini->paramnum );
	do newparnum = (long)GetValue(prompt, PARNUM_HE, (double)ini->paramnum);
	while( !Fail() && newparnum < 0 );
	if( Fail() ) return;
	if( newparnum == 0 )
	{
		ini->paramnum = 0;
		if( ini->param != NULL ) delete ini->param;
		ini->param = NULL;
		return;
	};

	// --- memory allocation -------------------------------
	if( ini->paramnum != newparnum || ini->param == NULL )
	{
		if( ini->param != NULL ) delete ini->param;
		ini->param = new double[ newparnum ];
		if( ini->param == NULL )
		{
			message( "No memory in configuration.",MSG_ERROR);
			FailFlag = 1==1;
			return;
		};
		ini->paramnum = newparnum;
		for( i=0; i < newparnum; ++i ) ini->param[i] = 0.0;
	};

	// --- reading values ----------------------------------
	for( i=0; !Fail() && i < newparnum; ++i )
	{
		sprintf( prompt, "%ld. constant parameter (%lf): ",
			i+1, ini->param[i] );
		ini->param[i] = GetValue( prompt, "", ini->param[i] );
	};
};


// -----------------------------------------------------------------------


void	Configure::GetBoundVectors( long newdim ) {

	long	i;
	double	d;

	if( Fail() ) return;
	if( !BOUND_NEEDED(ini->type) )
	{
		if( ini->lowb != NULL ) delete ini->lowb;
		if( ini->upb != NULL ) delete ini->upb;
		ini->upb = ini->lowb = NULL;
		return;
	};

	// --- memory allocation -------------------------------
	if( ini->dimension != newdim || ini->upb == NULL || ini->lowb == NULL )
	{
		if( ini->lowb != NULL ) delete ini->lowb;
		if( ini->upb != NULL ) delete ini->upb;
		ini->lowb = new double[ newdim ];
		ini->upb = new double[ newdim ];
		if( ini->upb == NULL || ini->lowb == NULL )
		{
			message( "No memory in configuration.",MSG_ERROR);
			FailFlag = 1==1;
			return;
		};
		for( i=0; i < newdim; ++i ) ini->upb[i] = ini->lowb[i] = 0.0;
	};

	// --- reading values ----------------------------------
	for( i=0; !Fail() && i < newdim; ++i )
	{
		sprintf( prompt, "%ld. component of lower bound (%lf): ",
			i+1, ini->lowb[i] );
		ini->lowb[i] = GetValue( prompt, "", ini->lowb[i] );
	};
	for( i=0; !Fail() && i < newdim; ++i )
	{
		sprintf( prompt, "%ld. component of upper bound (%lf): ",
			i+1, ini->upb[i] );
		do d = GetValue( prompt, "", ini->upb[i] );
		while( !Fail() && d < ini->lowb[i] );
		ini->upb[i] = d;
	};
};


// -----------------------------------------------------------------------
// --- Functions that do automatic parameter setting
// -----------------------------------------------------------------------


long	Configure::Diff() {
// this should be made 0 by the parameter setting algorithm

	unsigned long	evals1, evals2, i;

	evals1 = ini->evals[0];
	evals2 = ini->maxevals;
	for( i=1; i < ini->levels; ++i )
	{
		evals1 += ini->evals[i];
		evals2 -= ini->newspecevals[i-1];
	};

	return evals1 - evals2;
};


// -----------------------------------------------------------------------


void	Configure::SetMissingN() {
// everything is set when it is called

	if( Fail() ) return;

	ini->maxevals = ini->evals[0];
	for( long i = 1; i<ini->levels; ++i )
		ini->maxevals += ini->evals[i] + ini->newspecevals[i-1];
};


// -----------------------------------------------------------------------


void	Configure::SetMissingNu() {

	double		evals1; // double to handle negative numbers
	unsigned long	evals2, i;

	if( Fail() ) return;

	evals1 = ini->maxevals;
	for( i=1; i < ini->levels; ++i ) evals1 -= ini->newspecevals[i-1];
	if( evals1 < 0 )
	{
		message("All evaluations (N) too small, parameters not set.",
								MSG_ERROR);
		FailFlag = 1==1;
		return;
	};

	ini->threshold = 1.0;
	SetVectors( NOTHING_MISSING );
	evals2 = 0;
	for( i=0; i < ini->levels; ++i ) evals2 += ini->evals[i];
	ini->threshold = evals1 / (double)evals2;
	SetVectors( NOTHING_MISSING );
};


// -----------------------------------------------------------------------


void	Configure::SetMissingLastR() {
// the interval for last_r is (0,r[0]); Diff is strictly decreasing in last_r
// on this interval

	double	left, right, center;
	long	diff;

	if( Fail() ) return;

	// --- feasibility condition ---------------------
	last_r = ini->r[0];
	SetVectors( NOTHING_MISSING );
	if( Diff() > 0 )
	{
		message("No solution, parameters not set.",MSG_ERROR);
		FailFlag = 1==1;
		return;
	};

	// --- approximating last_r ------------------------
	left = 0;
	right = ini->r[0];
	for( long i = 0; i < 1000; ++i )
	{
		center = ( left + right ) / 2.0;
		last_r = center;
		SetVectors( NOTHING_MISSING );
		diff = Diff();
		if( diff < 0 ) right = center;
		else if( diff > 0 ) left = center;
		else break;
	};
};


// -----------------------------------------------------------------------


void	Configure::SetMissingL() {
// the interval for L is (2,+inf); Diff goes to inf as L goes to inf

	long	olddiff, diff;

	if( Fail() ) return;

	// --- feasibility condition ---------------------
	ini->levels = 2;
	ReAllocVectors( 2 );
	if( Fail() ) return;
	SetVectors( NOTHING_MISSING );
	if( (diff = Diff()) > 0 )
	{
		message("No solution, parameters not set.",MSG_ERROR);
		FailFlag = 1==1;
		return;
	};

	// --- approximating L ------------------------
	do
	{
		ReAllocVectors( ++ini->levels );
		if( Fail() ) return;
		SetVectors( NOTHING_MISSING );
		olddiff = diff;
		diff = Diff();
	}
	while( diff < 0 );
	if( -olddiff < diff )
	{
		ReAllocVectors( --ini->levels );
		if( Fail() ) return;
		SetVectors( NOTHING_MISSING );
	};
};


// -----------------------------------------------------------------------


void	Configure::SetMissingM() {
// the interval for M is (1,+inf); Diff is strictly increasing in M
// on this interval

	double	left, right, center;
	long	diff;

	if( Fail() ) return;

	// --- feasibility condition ---------------------
	ini->maxspecnum = 1;
	SetVectors( NOTHING_MISSING );
	if( Diff() > 0 )
	{
		message("No solution, parameters not set.",MSG_ERROR);
		FailFlag = 1==1;
		return;
	};

	// --- approximating M ------------------------
	left = 1;
	do
	{
		right = (ini->maxspecnum *= 2);
		SetVectors( NOTHING_MISSING );
	}
	while( Diff() < 0 );
	for( long i = 0; i < 1000; ++i )
	{
		center = ( left + right ) / 2.0;
		ini->maxspecnum  = (long) center;
		SetVectors( NOTHING_MISSING );
		diff = Diff();
		if( diff < 0 ) left = center;
		else if( diff > 0 ) right = center;
		else break;
	};
};


// -----------------------------------------------------------------------


void	Configure::SetVectors( long miss ) {

	double	beta;
	long	i;

	// --- new species evaluations -----------------------------------
	if( miss != MISSING_M )
		for( i=1; i < ini->levels; ++i )
			ini->newspecevals[i-1] = 3 * ini->maxspecnum;

	// --- radius vector ---------------------------------------------
	if( miss != MISSING_LAST_R )
	{
		// --- calculate exponential factor for radii -------------
		beta = pow( (last_r / ini->r[0]), 1.0/(ini->levels-1) );

		// --- calculate radii ------------------------------------
		for( i=1; i < ini->levels; ++i )
			ini->r[i] = ini->r[i-1] * beta;
	};

	// --- evaluations --------------------------------------------------
	if( miss != MISSING_M && miss != MISSING_NU && miss != MISSING_LAST_R )
		for( i=1; i < ini->levels; ++i )
			ini->evals[i] = (unsigned long) (
				Master::M( i, ini ) * ini->r[0] *
				ini->threshold /
				ini->prototype->v( ini->r[i] ) );
};


// -----------------------------------------------------------------------


void	Configure::ReAllocVectors( long levels ) {

	FailFlag = 1==1;

	// --- memory allocation --------------------------------------
	if( ini->evals != NULL ) delete ini->evals;
	if( ini->newspecevals != NULL ) delete ini->newspecevals;
	if( ini->r != NULL ) delete ini->r;
	ini->evals = new unsigned long[ levels ];
	ini->newspecevals = new unsigned long[ levels - 1 ];
	ini->r = new double[ levels ];
	if( ini->evals==NULL || ini->newspecevals==NULL || ini->r==NULL )
	{
		message("No memory for storing parameters.",MSG_ERROR);
		return;
	};

	// --- method independent defaults ----------------------------
	ini->r[0] = ini->prototype->Diameter( ini ); // diameter of space
	ini->evals[0] = 0;

	FailFlag = 1==0;
};


// -----------------------------------------------------------------------


void	Configure::AutoParameters( long miss ) {

	if( Fail() ) return;

	// --- initializations ----------------------------------------
	if( miss == MISSING_L ) ini->levels = 2;
	if( ini->levels == 1 )
	{
		ReAllocVectors( 2 );
		if( Fail() ) return;
		ini->evals[0] = ini->maxevals;
		goto WARN;
	};
	ReAllocVectors( ini->levels );
	if( Fail() ) return;
	SetVectors( miss ); // sets everything possible

	// --- calculating missing parameter --------------------------
	switch( miss )
	{
		case MISSING_N:
			SetMissingN();
			break;
		case MISSING_NU:
			SetMissingNu();
			break;
		case MISSING_LAST_R:
			SetMissingLastR();
			break;
		case MISSING_L:
			SetMissingL();
			SetMissingN();
			break;
		case MISSING_M:
			SetMissingM();
			break;
	};
	
	// --- send warnings if needed ----------------------------------
WARN:	if( ini->levels == 1 )
	{
		message( "levels set to 1, special uego things won't be used!",
							MSG_INFORMATION);
		return;
	};
	if( ini->maxspecnum > Master::M( ini->levels-1, ini ) )
	{
		sprintf( prompt, "Effective max. species number is only %ld!",
				(long) Master::M( ini->levels-1, ini ) );
		message( prompt, MSG_INFORMATION );
	};
	for( long i=1; i<ini->levels; ++i )
		if( ini->Evals(i) / (long)Master::M( i, ini ) == 0 )
		{
			sprintf( prompt,
			"0 evaluations for level %ld species!", i+1 );
			message( prompt, MSG_INFORMATION );
		};
};


// -----------------------------------------------------------------------
// --- Functions for configuration
// -----------------------------------------------------------------------


void	Configure::DisplayConfig() {

	long	l;
	double	d;

	FailFlag = 1==0;
	printf( STARTUP_HE );

	// --- random numbers ---------------------------------------------
	sprintf( prompt, "%s (%lu): ", SEED_T, ini->seed );
	ini->seed = (unsigned long)GetValue(prompt, SEED_HE, (double)ini->seed);

	// --- objective function ------------------------------------------
	sprintf( prompt, "%s (%ld): ", TYPE_T, ini->type );
	do l = (long) GetValue( prompt, TYPE_HE, (double)ini->type );
	while( !Fail() && !VALID_TYPE(l) );
	ini->type = l;

	sprintf( prompt, "%s (%ld): ", FNUM_T, ini->fnum );
	ini->fnum  = (long) GetValue( prompt, FNUM_HE, (double)ini->fnum );

	sprintf( prompt, "%s (%ld): ", DIM_T, ini->dimension );
	do l = (long) GetValue( prompt, DIM_HE, (double)ini->dimension );
	while( !Fail() && l < 1 );

	GetBoundVectors( l );
	ini->dimension = l;

	GetParameterVector(); // reads also paramnum

	ini->SetPrototype();
	FailFlag == FailFlag || ini->Fail();

	// --- uego parameters ---------------------------------------------
	sprintf( prompt, "%s (%lu): ", MAXEV_T, ini->maxevals );
	ini->maxevals = (unsigned long) GetValue( prompt, MAXEV_HE,
							(double)ini->maxevals );

	sprintf( prompt, "%s (%ld): ", MAXSPEC_T, ini->maxspecnum );
	do l = (long) GetValue( prompt, MAXSPEC_HE, (double)ini->maxspecnum );
	while( !Fail() && l < 1 );
	ini->maxspecnum = l;

	sprintf( prompt, "%s (%lg): ", THR_T, ini->threshold );
	do d = GetValue( prompt, THR_HE, ini->threshold );
	while( !Fail() && d <= 0.0 );
	ini->threshold = d;

	sprintf( prompt, "%s (%ld): ", LEVELS_T, (long)ini->levels );
	do l = (long) GetValue( prompt, LEVELS_HE, (double)ini->levels );
	while( !Fail() && l < 1 );
	ini->levels = (short)l;

	if( !Fail() && ini->levels == 1 )
	{
		AutoParameters( NOTHING_MISSING );
		return;
	};

	// --- automatic parameters ---------------------------------------
	sprintf( prompt, "%s (%lg): ", LAST_R_T, last_r );
	do d = GetValue( prompt, LAST_R_HE, last_r );
	while( !Fail() && (d <= 0.0 || d >= ini->prototype->Diameter(ini)));
	last_r = d;

	sprintf( prompt, "Parameter to ignore (2): " );
	do l = (long) GetValue( prompt, MISS_HE, 2.0 );
	while( !Fail() && !VALID_MISS(l) );

	if( !Fail() ) AutoParameters( l ); // sets the remaining vectors
};


// -----------------------------------------------------------------------


void	Configure::ParamConfig( int argc, char** argv ) {
// the parameter order is fixed in command line

	unsigned long	ul;
	long		l, i=0;
	double		d;
	FILE*		parfile;

	FailFlag = 1==0;

	// --- random numbers ---------------------------------------------
	if( !FailFlag && i < argc && argv[i][1]=='s' )
	{
		if( sscanf( argv[i++]+2, "%lu", &ul ) == 1 ) ini->seed = ul;
		else FailFlag = 1==1;
	};

	// --- objective function ------------------------------------------
	if( !FailFlag && i < argc && argv[i][1]=='p' )
	{
		parfile = fopen( argv[i++]+2, "rt" );
		if( parfile != NULL )
		{
			ReadParFile( parfile );
			fclose( parfile );
		}
		else FailFlag = 1==1;
	};

	ini->SetPrototype();
	FailFlag == FailFlag || ini->Fail();

	// --- uego parameters ---------------------------------------------
	if( !FailFlag && i < argc && argv[i][1]=='N' )
	{
		if( sscanf( argv[i++]+2, "%lu", &ul ) == 1 ) ini->maxevals = ul;
		else FailFlag = 1==1;
	};

	if( !FailFlag && i < argc && argv[i][1]=='M' )
	{
		if( sscanf( argv[i++]+2, "%ld", &l ) == 1 && l >= 1 )
			ini->maxspecnum = l;
		else FailFlag = 1==1;
	};

	if( !FailFlag && i < argc && argv[i][1]=='t' )
	{
		if( sscanf( argv[i++]+2, "%lf", &d ) == 1 && d >= 0.0 )
			ini->threshold = d;
		else FailFlag = 1==1;
	};

	if( !FailFlag && i < argc && argv[i][1]=='l' )
	{
		if( sscanf( argv[i++]+2, "%ld", &l ) == 1 && l >= 1 )
			ini->levels = (short)l;
		else FailFlag = 1==1;
	};

	// --- automatic parameters ---------------------------------------
	if( !FailFlag && i < argc && argv[i][1]=='r' )
	{
		if( sscanf( argv[i++]+2, "%lf", &d ) == 1 &&
			d > 0.0 && d < ini->prototype->Diameter(ini) )
			last_r = d;
		else FailFlag = 1==1;
	};

	l = 2;
	if( !FailFlag && i < argc && argv[i][1]=='x' )
		if( sscanf( argv[i++]+2, "%ld", &l ) != 1 || !VALID_MISS(l) )
			FailFlag = 1==1;

	if( !FailFlag ) AutoParameters( l ); // sets the remaining vectors
	else message("Command line configuration failed.",MSG_ERROR);
};


// -----------------------------------------------------------------------


void	Configure::Config( char* ininame, int argc, char** argv ) {

	FILE	*inifile = NULL;

	FailFlag = 1==1;

	// --- creating default Ini object --------------------------------
	inifile = fopen( ininame, "rt" );
	if( inifile == NULL )
		ini = new Ini;	// default ini
	else
		ini = new Ini( inifile );
	if( ini == NULL || ini->Fail() )
	{
		message("Could not create default settings.",MSG_ERROR);
		return;
	};
	last_r = ini->R( ini->Levels()-1 );

	// --- configuring ------------------------------------------------
	if( argc > 0 ) ParamConfig( argc, argv );
	else DisplayConfig();

	// --- exiting ----------------------------------------------------
	if( inifile != NULL ) fclose( inifile );
	if( !Fail() )
	{
		inifile = fopen( ininame, "wt");
		ini->Save( inifile );
		FailFlag = ini->Fail();
	};
	if( !Fail() )
		message( "Configuration successful.", MSG_INFORMATION );
	delete ini;
};
