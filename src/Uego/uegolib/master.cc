#include <string.h>
#include "uego.h"


// -----------------------------------------------------------------------


void	Master::ReInit( Ini* ini, char* trace, char *file) {

	SearchSpElement	*root;

	FailFlag = 1==1;

	

	_ini = ini;		// only storing pointer!
	tracename = trace;	// only storing pointer!

	
	
	// --- initializing species list with root:
	// --- 1 random element from space
	head = new SpeciesList;
	if( head != NULL ) root = _ini->Prototype()->RandNew();
	if( root != NULL ) head->next = new SpeciesList( root, 0 );
	if( head == NULL || root == NULL || head->next == NULL )
	{
		message((char*)"No memory for master.",MSG_ERROR);
		return;
	};

	BestObj = root->CurrValue();
	BestX = head->next;
	
	head->next->prev = head;
	length = 1;
	FunctionEvals = 1;
	
	//printf(" %e %ld \n", BestObj, FunctionEvals); 
	SaveBest(file);

	
	level = 0;
	FailFlag = 1==0;
};


// -----------------------------------------------------------------------


void	Master::Clean() {

	SpeciesList	*tmp;

	if( head != NULL )
	{
		while( head->next != NULL )
		{
			tmp = head->next;
			head->next = head->next->next;
			delete tmp;
		};
		delete head;
	};

	// _ini and tracename must be deleted by the caller of the constructor
};


// -----------------------------------------------------------------------


void	Master::Fuse() {

	const double	r = _ini->R( level );
	//printf("Level %i radio %lf \n", level, r);
	//getchar();
	SpeciesList	*tmp, *act;

	if( length < 2 ) return;

	for( act = head->next->next; act != NULL; act = act->next )
	{
		tmp = act->prev;
		while( tmp != head )
		{

			// --- compare act and tmp for fusion
			if( act->center->Distance( tmp->center ) < r )
			{
				if( act->level <= tmp->level )
				{
					// --- absorbing tmp
					act->Absorb( tmp );

					// --- cutting out tmp
					// (tmp is never the last!)
					length--;
					tmp = tmp->prev;
					tmp->next = tmp->next->next;
					delete tmp->next->prev;
					tmp->next->prev = tmp;
				}
				else
				{
					// --- absorbing act
					tmp->Absorb( act );

					// --- cutting out act
					length--;
					tmp = act; // old tmp not used anymore
					act = act->prev;
					act->next = act->next->next;
					delete tmp;
					if( act->next != NULL )
						act->next->prev = act;
					break;
				};
			}
			else tmp = tmp->prev; // nothing has been cut out
		};
	};
};


// -----------------------------------------------------------------------


void	Master::_NewSpecies( long evals, char *file ) {

	SpeciesList	newspechead,	// First element is newspechead.next!
			*newlst,
			*newend,	// shows end of new spec list
			*tmp;		// position in old list
	long 		newevals;

	FailFlag = 1==1;

	newend = &newspechead;
	for( tmp = head; tmp->next != NULL; tmp = tmp->next )
	{
		// --- to allow multiple calls w.o. level increased
		if( tmp->next->level == level ) continue;


		// --- list of new species in newlst
		newevals = tmp->next->NewSpecies( &newlst, level, evals );
		FunctionEvals +=newevals;
		if( tmp->next->Fail() ) return;
		
		//---Check the best
		//---El centro se compara con todos los puntos, con lo que siempre tiene el mejor.
		if( tmp->next->center->CurrValue() > BestObj ) {
			BestObj = tmp->next->center->CurrValue();
			BestX = tmp->next;
			SaveBest(file);
		}


		// --- insert new list to newspechead, increase length
		newend->next = newlst;
		if( newlst != NULL ) newlst->prev = newend;
		while( newend->next != NULL )
		{
			++length;
			newend = newend->next;
			
		};
	};
	tmp->next = newspechead.next;
	if( newspechead.next != NULL ) newspechead.next->prev = tmp;
	newspechead.next = NULL; // to prevent destructing;

	FailFlag = 1==0;
};


// -----------------------------------------------------------------------


void	Master::NewSpecies(char *file) {

	const long	oldlength = length;
	SpeciesList *tmp;
	
	_NewSpecies( _ini->NewSpecEvals(level) / length ,file);
	
	// ---------------------------------------------------
	// --- beginning of species creation forcing extension
	// ---------------------------------------------------
	if( Fail() ) return;
	
	/*const long	needed_msn = (long)pow( M( _ini->Levels() - 1 ),
					(double)level / (_ini->Levels()-1)),
			maxevals = head->Evals(level) / oldlength; 
	long		max_spec_num = (long)M( level );	
	*/
	
	long	max_spec_num = (long)M( level );	
	long	bugetSpecies = INI.Evals(level-1) / max_spec_num ; 		//---Presupuesto que tenía cada especie en el nivel anterior para optimizar
	long	maxevals = INI.Evals(level-1)-(bugetSpecies*oldlength); 	//---Evaluaciones remanentes

	/*printf("Nivel actual = %i \n", level);
	printf("Longitud inicial = %ld \n", oldlength);
	printf("Evals en el nivel anterior %ld",  INI.Evals(level-1)); 
	printf("Presupuesto por especie en el nivel %i %ld \n", level-1, bugetSpecies);
	printf("Evaluaciones remanentes %ld \n", maxevals);
	getchar();*/


	Fuse();	if( Fail() || maxevals < 3 ) return;
	CheckLength();
	
	//if( maxevals < 3 ) return;
	//printf("Estoy aquí\n");
	_NewSpecies( maxevals, file); if( Fail() ) return;  
	
	//Fuse(level);		 if( Fail() ) return;  
	//--max_spec_num;
	//CheckLength( max_spec_num );
	
	// ---------------------------------------------------
	// --- end of species creation forcing extension
	// ---------------------------------------------------
};


// -----------------------------------------------------------------------


double	Master::M( long i, Ini* ini ) {

	return ini->MaxSpecNumber();
};

// -----------------------------------------------------------------------


void	Master::Optimize(char *file) {

	SpeciesList	*tmp;
	long		newevals;
	
	long	max_spec_num = (long)M( level );	
	long	bugetSpecies = INI.Evals(level) / max_spec_num ; //---Presupuesto por especie en el nivel actual
	

	FailFlag = 1==1;
					
	for( tmp = head->next; tmp != NULL; tmp = tmp->next )
	{
		
		newevals=tmp->Optimize(bugetSpecies);
		//printf("newevals= %ld \n", newevals);
		//getchar();
		FunctionEvals+=newevals;

		if( tmp->center->CurrValue() > BestObj ) {
			BestObj = tmp->center->CurrValue();
			BestX = tmp;
			SaveBest(file);
		}

		if( tmp->Fail() ) return;
	};

	FailFlag = 1==0;
};


// -----------------------------------------------------------------------


void	Master::CheckLength( long new_length ) {

	SpeciesList	*end;

	if( new_length == -1 ) new_length = _ini->MaxSpecNumber();
	if( length > new_length )
	{
		message((char*) "Too many species, shortening species list.",
							MSG_INFORMATION);
		end = head;
		while( end->next != NULL ) end = end->next;

		// --- deleting from end of list
		for(; length > new_length; --length )
		{
			end = end->prev;
			delete end->next;
		};
		end->next = NULL;
	};
};


// -----------------------------------------------------------------------


void	Master::_Go(char * file) {

	//FunctionEvals += head->next->Optimize( _ini->Evals(0) );

	SaveGeneration(file);

	SpeciesList *tmp;
	for( level=1; level < _ini->Levels(); ++level )
	{
		
		//printf("\n \n LEVEL = %i \n \n", level);
		//---Creation
		/*printf("%----------------------------------------------------------------%\n");
		printf("%-------------------------NewSpecies-----------------------------%\n");
		printf("%----------------------------------------------------------------%\n");*/
		NewSpecies(file);	if( Fail() ) return;
		
		
		/*printf("\n Longitud lista= %ld::\n",length);
		for( tmp = head; tmp->next != NULL; tmp = tmp->next )
		{
			double x[30];
			tmp->next->center->GetX(x);
			printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  				Obj= %lf %ld\n",x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],
			x[20],x[21],x[22],x[23],x[24],x[25],x[26],x[27],x[28],x[29],tmp->next->center->CurrValue(),tmp->next->level);
		}*/
		
		//---Fusion
		/*printf("%----------------------------------------------------------------%\n");
		printf("%-------------------------Fusion---------------------------------%\n");
		printf("%----------------------------------------------------------------%\n");
		*/
		Fuse();	if( Fail() ) return;
		
		/*printf("\n Longitud lista= %ld::\n",length);
		for( tmp = head; tmp->next != NULL; tmp = tmp->next )
		{
			double x[30];
			tmp->next->center->GetX(x);
			printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  				Obj= %lf %ld\n",x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],
			x[20],x[21],x[22],x[23],x[24],x[25],x[26],x[27],x[28],x[29],tmp->next->center->CurrValue(),tmp->next->level);
		}*/


		//---CheckLength
		/*printf("%----------------------------------------------------------------%\n");
		printf("%-------------------------CheckLength----------------------------%\n");
		printf("%----------------------------------------------------------------%\n");*/

		CheckLength();

		/*printf("\n Longitud lista= %ld::\n",length);
		for( tmp = head; tmp->next != NULL; tmp = tmp->next )
		{
			double x[30];
			tmp->next->center->GetX(x);
			printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  				Obj= %lf %ld\n",x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],
			x[20],x[21],x[22],x[23],x[24],x[25],x[26],x[27],x[28],x[29],tmp->next->center->CurrValue(),tmp->next->level);
		}*/
		

		//---Optimization

		/*printf("%----------------------------------------------------------------%\n");
		printf("%-------------------------Optimization---------------------------%\n");
		printf("%----------------------------------------------------------------%\n");*/

		Optimize(file);	if( Fail() ) return;

		/*printf("\n Longitud lista= %ld::\n",length);
		for( tmp = head; tmp->next != NULL; tmp = tmp->next )
		{
			double x[30];
			tmp->next->center->GetX(x);
			printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  				Obj= %lf %ld\n",x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],
			x[20],x[21],x[22],x[23],x[24],x[25],x[26],x[27],x[28],x[29],tmp->next->center->CurrValue(),tmp->next->level);
		}*/
		
		//---Fusion
		Fuse();		if( Fail() ) return;

		SaveGeneration(file);
	};
};

// -----------------------------------------------------------------------

void	Master::SaveBest( char * file) {

	
	FailFlag = 1==1;

	
	
	
	FILE *fptable1;
	char tablefile1[90];
	strcpy(tablefile1,"BESTOBJ_VS_FUNCEVALU");	  
	strcat(tablefile1,file);


	FILE *fptable2;
	char tablefile2[90];
	strcpy(tablefile2,"BEST-SOL_");	  
	strcat(tablefile2,file);

	
	//---Saving best objective function value in each level (to depict evolution along the generations)
	if ((fptable1 = fopen(tablefile1,"a+")) == 0) printf("I can't open/create the file 1: %s\n",tablefile1);
	fprintf( fptable1, "%.10le   %ld\n ", BestObj, FunctionEvals );

	//---Solution associated to the best objective function value. 	
	if ((fptable2 = fopen(tablefile2,"a+")) == 0) printf("I can't open/create the file 1: %s\n",tablefile2);
	BestX->Save2( fptable2 );

	fclose(fptable1);
	fclose(fptable2);

	FailFlag = 1==0;
};

// -----------------------------------------------------------------------

void	Master::SaveGeneration( char * file) {

	SpeciesList	*tmp,*tmp_best;
	double		best, best_x0,best_x1;
	
	FailFlag = 1==1;

	
	// --- saving pareto points
	
	FILE *fptable1;
	char tablefile1[90];
	strcpy(tablefile1,"EVOLUTION-BEST-OBJ_PER_ITER");	  // --- Poblaciones en cada nivel
	strcat(tablefile1,file);


	FILE *fptable2;
	char tablefile2[90];
	strcpy(tablefile2,"POPULUATIONS_PER_ITER");	  // --- Poblaciones en cada nivel
	strcat(tablefile2,file);


	FILE *fptable3;
	char tablefile3[90];
	strcpy(tablefile3,"BEST-SOL_PER_ITER"); 	   //--- Mejor solución encontrada en cada nivel
	strcat(tablefile3,file);


	FILE *fptable4;
	char tablefile4[90];
	strcpy(tablefile4,"OBJ_FUNC_PER_ITER"); 	   //---Número de evaluaciones consumidas por nivel
	strcat(tablefile4,file);



	tmp = head->next;
	tmp_best=head->next;
	best = tmp->center->CurrValue();
	for( ; tmp != NULL; tmp = tmp->next )
		if( tmp->center->CurrValue() > best ){
			best = tmp->center->CurrValue();
			tmp_best=tmp;
		   }

	//---Saving best objective function value in each level (to depict evolution along the generations)
	if ((fptable1 = fopen(tablefile1,"a+")) == 0) printf("I can't open/create the file 1: %s\n",tablefile1);
	fprintf( fptable1, "%.10le\n ", best );


	// --- Saving best individual found in each level. For each one of them, the active RBFs, the centers and the sigmas are showed. 
	// --- Additionally, the number of active RBFs, the minimum, average, maximum and standard deviation of all the signmas values is computed

	if ((fptable3 = fopen(tablefile3,"a+")) == 0) printf("I can't open/create the file 1: %s\n",tablefile3);
	
	
	fprintf( fptable3, " -------------------------------- \n " );
	fprintf( fptable3, "---------- Level = %i ---------- \n ", level );
	fprintf( fptable3, "-------------------------------- \n " );

	tmp_best->Save2(fptable3);	
	if( tmp_best->Fail() ) return;
	

	
	// --- Saving the whole population in each level. 
	// --- For each individual in the population the active RBFs, the centers and the sigmas are showed. 
	// --- Additionally, the number of active RBFs, the minimum, average, maximum and standard deviation of all the signmas values is computed

	if ((fptable2 = fopen(tablefile2,"a+")) == 0) printf("I can't open/create the file 1: %s\n",tablefile2);
	
	fprintf( fptable2, "-------------------------------- \n " );
	fprintf( fptable2, "---------- Level = %i ---------- \n ", level );
	fprintf( fptable2, "-------------------------------- \n " );

	tmp = head->next;
	while( tmp != NULL )
	{       tmp->Save2( fptable2 );
		if( tmp->Fail() ) return;
		tmp = tmp->next;
	};
	fprintf( fptable2, " \n \n " );
	

	//---Saving the number of function evaluation on each level
	if ((fptable4 = fopen(tablefile4,"a+")) == 0) printf("I can't open/create the file 1: %s\n",tablefile4);
	fprintf( fptable4, "%ld\n ", FunctionEvals );



	fclose(fptable1);
	fclose(fptable2);
	fclose(fptable3);
	fclose(fptable4);


	FailFlag = 1==0;
};


// -----------------------------------------------------------------------


void	Master::Save( FILE* stream, double tiempo,unsigned long flags ) {

	SpeciesList	*tmp,*tmp_best;
	double		best, best_x0,best_x1;

	FailFlag = 1==1;

	if( flags & SAVE_INI )
	{
		_ini->Save( stream );
		if( _ini->Fail() ) return;
	};

	if( fprintf( stream,
		"%ld %ld %lf ",FunctionEvals, length,tiempo ) == EOF )
	{
		message((char*)"Error saving master.",MSG_ERROR);
		return;
	};

	// --- which is best?
	tmp = head->next;
	tmp_best=head->next;
	best = tmp->center->CurrValue();
	for( ; tmp != NULL; tmp = tmp->next )
		if( tmp->center->CurrValue() > best ){
			best = tmp->center->CurrValue();
			tmp_best=tmp;
		   }
	tmp_best->Save( stream );	
	if( tmp_best->Fail() ) return;
		
	// --- saving species
	printf("\n \n Total Species List with conversion in first elements, i.e. if x>0.5 then x=1\n");	
	if( !(flags & SHORT_SAVE) )
	{
		tmp = head->next;
		while( tmp != NULL )
		{       tmp->Save( stream );
			printf("\n");
			if( tmp->Fail() ) return;
			tmp = tmp->next;
		};
	};

	FailFlag = 1==0;
};

// -----------------------------------------------------------------------
void	Master::sortObj() {
	SearchSpElement *root;
	short auxLevel;
	SpeciesList	*tmp;
	root = _ini->Prototype()->RandNew();
	
	for( long i = 0; i< length; i++) 
	{   	
	    tmp = head->next;
	    while(tmp->next!=NULL)
	    {
		if( tmp->center->CurrValue() < tmp->next->center->CurrValue() ) 
		{
			// --- sort
			root->UpdateFrom(tmp->center);
			auxLevel = tmp->level;

			tmp->center->UpdateFrom(tmp->next->center);
			tmp->level = tmp->next->level;
			
			tmp->next->center->UpdateFrom(root);
			tmp->next->level = auxLevel;			

		}//end if
		tmp = tmp->next;
	   };//end while
	};//end for

}


// -----------------------------------------------------------------------


void	Master::CheckPoint() {

	char	name[500], msg[500];
	FILE*	outf;

	FailFlag = 1==0;
	if( tracename == NULL ) return;
	else FailFlag = 1==1;

	strcpy( name, tracename );
	sprintf( name+strlen(name), "%04ld.ckp", tracenum );
	if( (outf = fopen( name, "wt" )) == NULL )
	{
		sprintf( msg, "Could not open file '%s'.", name );
		message( msg, MSG_ERROR );
		return;
	};

	Save( outf ); // Sets FailFlag
	fclose( outf );

	++tracenum;
};

