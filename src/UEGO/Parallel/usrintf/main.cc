#include "time.h"
#include <string.h>
#include "usrintf.h"
#include "uego.h"
#include "configur.h"



Ini*	Master::_ini = NULL; 	// initialization of static member
void	setMsgLevel( char );

//--- The following is a non-documented feature, to use it in a modul -----
//--- declare global variables as 'extern' --------------------------------

char**	RemArgv = NULL;	// command line pars after '--'
int	RemArgc = 0;	// num. of command line pars after '--' 


void	CutCommandLine( int* argc, char** argv ) {

	for( int i=0; i<*argc; ++i )
	{
		if( strcmp( argv[i], "--" ) == 0 )
		{
			RemArgc = *argc-(i+1);
			if( i != *argc-1 ) RemArgv = argv + (i+1);
			*argc = i;
			break;
		};
	};
};
		

//--- These are local to this modul ------------------------------------

unsigned long	saveflags;	// GetPars() sets it
long		repcount;	// GetPars() sets it
char		*trace,		// GetPars() sets it
		msg[250];	// working


//-----------------------------------------------------------------------


char	GetPars( int argc, char** argv ) {
// sets repcount and trace; returns false on failure

	trace = NULL;
	repcount = 1;
	saveflags = 0;

	for(int i=0; i<argc; ++i)
	{
		switch( argv[i][1] )
		{
		case 'T':
			if( trace != NULL ) delete trace;
			trace = new char[ strlen(argv[i]+2) + 1 ];
			if( trace == NULL )
			{
				message((char*)"No memory for comm. line.",MSG_ERROR);
				return 1==0;
			};
			strcpy( trace, argv[i]+2 );
			break;
		case 'r':
			if( sscanf( argv[i]+2, "%ld", &repcount ) != 1 )
			{
				message((char*)"Bad repeate count.",MSG_ERROR);
				return 1==0;
			};
			break;
		case 's':
			saveflags |= SHORT_SAVE;
			break;
		};
	};

	return 1==1;
};


//-----------------------------------------------------------------------


int	main( int argc, char** argv ) {

	int	errcode = 0;
	Ini	*ini = NULL;
	Master	*master = NULL;
	FILE	*inifile = NULL;
	double  tiemporeloj, time1, time2;
	time_t t1,t2;
	clock_t clock1,clock2;

	//-- Parallel variable
	int	namelen=0,
		myid=0, numproc=0;			
	char	name[25];
	

	//---------------------		
	MPI_Init(&argc,&argv);	
	MPI_Get_processor_name(name, &namelen);
	if(MPI_Comm_rank(MPI_COMM_WORLD,&myid)!=MPI_SUCCESS) 		
		printf("An error in: Ini::INI --> MPI_Comm_Rank\n");
	
	if (MPI_Comm_size(MPI_COMM_WORLD,&numproc)!=MPI_SUCCESS)	
		printf("An error in: Ini::INI --> MPI_Comm_size\n");
	//---------------------		
	


	setMsgLevel( MSG_INFORMATION );

	CutCommandLine( &argc, argv );	// cuts things after '--'

	if( argc < 2 )
	{
		sprintf( msg, "Version: %s                              Contact: jelasity@usa.net\n\n",UEGO_VERSION);
		message( msg, MSG_NOTHING );
		message((char*) DEFAULT_MESSAGE, MSG_NOTHING );
		return 0;
	}
	else if( argv[1][0] != '-' ) // assuming ini filename extension
	{
		sprintf( msg, "Opening ini file %s.%s",UEGO_ININAME,argv[1]);
		message( msg, MSG_INFORMATION );
		sprintf( msg, "%s.%s", UEGO_ININAME, argv[1] );
		inifile = fopen( msg, "rt" );
	}
	else if( argv[1][1] != 'c' )
	{
		message( (char*)"Bad or unknown command line parameters.",MSG_ERROR);
		return 6;
	}
	else  // configure
	{
		sprintf( msg, "%s.%s", UEGO_ININAME, argv[1]+2 );
		Configure::Config( msg, argc-2, argv+2 );
		if( Configure::Fail() ) return 7;
		else return 0;
	};

	// --- Starting optimalization -------------------------------------
	message((char*)"Reading ini file.",MSG_INFORMATION);
	ini = new Ini( inifile );
	if( ini == NULL || ini->Fail() )
	{
		message((char*)"Could not read ini file.",MSG_ERROR);
		errcode = 1;
	}
	else
	{
		message((char*)"Creating master.",MSG_INFORMATION);
		if( GetPars( argc-2, argv+2 ) )
			master = new Master( ini, trace,argv[1]);
		
		if( master == NULL || master->Fail() )
		{
			message((char*)"Could not create master.",MSG_ERROR);
			errcode = 2;
		}
		else{
			message((char*)"Starting optimalization.",MSG_INFORMATION);
			// --- doing repcount experiments -----------------------
			for( long i=1; i<=repcount; ++i )
			{
				time1 = MPI_Wtime();
				t1=time(NULL);
				
				//clock1=clock();
				sprintf( msg, "Starting experiment %ld.", i );
				message( msg,MSG_INFORMATION );
				master->Go(argv[1]);
				time2 = MPI_Wtime();
				t2=time(NULL);
				
				//clock2=clock();
				tiemporeloj=(double)t2-t1;
				//tiemporeloj=(double)(clock2-clock1)/(double)CLOCKS_PER_SEC; 
				
				if( master->Fail() )
				{
					errcode = 3;
					break;
				}
				else
				{
					message((char*)"Saving results.",MSG_INFORMATION);
					master->Save( stdout,tiemporeloj,saveflags );
					if( master->Fail() )
					{
						errcode = 4;
						break;
					};
				};
			};
		}
	};

	if( errcode == 0 ) message((char*)"Success.",MSG_INFORMATION);

	if( master != NULL ) delete master;
	if( ini != NULL ) delete ini;
	if( trace != NULL ) delete trace;
	

	MPI_Barrier( MPI_COMM_WORLD );	
	MPI_Finalize();

	return errcode;
};
