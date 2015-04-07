#ifndef MASTER_H
#define MASTER_H

#include "wqueue.h"

class CommunicationManager;


////////////////////////////////////////////////////////////
// $Id: master.h,v 2.6 1998/03/29 10:39:46 jelasity Exp $
// master.h
// declares class Master;
////////////////////////////////////////////////////////////
// modification history:
//	Jelasity 98 01 24 new save flag
////////////////////////////////////////////////////////////

// --- flags for save --------------------------------------------------

#define SAVE_INI 1UL
#define SHORT_SAVE 2UL

// ---------------------------------------------------------------------

class Master {

private:

	static Ini	*_ini;		// parameters given by the user

	SpeciesList	*head;		// first element is head->next!
	long		length;		// length of species list
	long		_length;
	long		length_;
	
	long		FunctionEvals;	// number of function evaluations
	double		BestObj;	// best objective function value	
	SpeciesList	*BestX;		// Best solution
	
	short		level;		// actual strict level

	// Thread synchronization
	wqueue<SearchSpElement *> to_simulate_queue;
	wqueue<SearchSpElement *> simulated_queue;
	CommunicationManager* ComManager;

	bool 	end_simulation;
	bool	generated_species;

	void	CheckLength( long=-1 );	// if list too long, shortens it
	void	Fuse();			// fuses species list using 'level'
	void	_NewSpecies( long, char*);
	void	NewSpecies(char *);

	void	_NewSpeciesOptiParamCereb( long, char*);
	void	NewSpeciesOptiParamCereb(char *);
	void	Optimize(char *);

	// --- to trace optimization process
	long	tracenum;	// counts checkpoints
	char	*tracename;
	void	CheckPoint();

	// --- this is needed for performing more than 1 experiments -------
	void	ReInit( Ini*, char*,char*);
	void	ReInitParal( Ini*, char*, char *);
	void	Clean();
	void	NewSearch(char * file) { Clean();ReInitParal( _ini, tracename,file);}
	void	_Go(char * file); // performs the optimization process

	char	FailFlag;

public:

	Master( Ini* ini, char* trace,char *file ) {tracenum=0; ComManager=0; _ini = ini; tracename = trace; head = NULL; FailFlag = 1==0;};
	~Master() { Clean(); };

	void	Go(char * file) { NewSearch(_ini->LoadStateFile());  if( !Fail() ) _Go(file); };
	void	Save( FILE*,double= 0UL, unsigned long = 0UL ); // saves status of search
	void    SaveGeneration(char * file);
	void	SaveBest( char * file);
	void	SaveState(char * FileName); // Save the current status of the search


	static Ini&	ini() { return *_ini; }; // safe access to _ini
	static char	iniSet() { return _ini != NULL; };
	static double	M( long, Ini* = _ini );
		// max theoretically possible spec list length at given level

	char	Fail() { return FailFlag; };
	void	sortObj();


	//---Parallel methods
	
	void OptimizeParal(char *) ;
	void NewSpeciesParal(char *);
	void packingSingleSpecies(int, SpeciesList *,int);
	void sendList();
	void FinalizeParal();
	

};


#endif
