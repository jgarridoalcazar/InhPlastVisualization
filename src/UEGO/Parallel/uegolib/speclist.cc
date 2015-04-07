#include "uego.h"

#include <fstream>


// -----------------------------------------------------------------------

long	SpeciesList::Optimize( long maxevals )   { //short tmplevel,

	long	l;
	//level=tmplevel;
	if( maxevals == -1 )
	{

		maxevals = Evals(level);
		if( Fail() ) return 0;
	};

	/*double x[1],x2[1];
	center->GetX(x);
	printf("Antes   --> OPTIMIZE:: %lf Obj= %lf Maxevals= %ld, level= %ld\n",x[0],center->CurrValue(),maxevals,level);
	*/	
	

	l = center->Optimize( level, maxevals );
	
	/*center->GetX(x2);
	printf("Despues --> OPTIMIZE:: %lf Obj= %lf Maxevals= %ld, level= %ld\n",x2[0],center->CurrValue(),l,level);
	*/	
	FailFlag = center->Fail();
	return l;
};

// -----------------------------------------------------------------------

SearchSpElement	* SpeciesList::InitializeOptimize()   { //short tmplevel,


	return center->InitOptimizeParal(level);
};


// -----------------------------------------------------------------------

SearchSpElement *	SpeciesList::ResumeOptimize(SearchSpElement * tmp)   { //short tmplevel,


	return center->ResumeOptimize( tmp );
};



// -----------------------------------------------------------------------


long	SpeciesList::Evals( long _level ) {

	FailFlag = 1==1;

	if( !Master::iniSet() )
	{
		message((char*)"Ini not set in SpeciesList::Evals()",MSG_ERROR);
		return 0;
	};

	if( _level == -1 ) _level = level; // called without param
	if( _level == 0 ) _level = 1;
	
	FailFlag = 1==0;


	//printf("Master::ini().Evals(_level)= %ld (long)Master::M(_level) %ld ",Master::ini().Evals(_level),(long)Master::M(_level)); 
	//printf("evals = %ld \n", Master::ini().Evals(_level) / (long)Master::M(_level));
	//getchar();
	return  Master::ini().Evals(_level) / (long)Master::M(_level);

	
};


// -----------------------------------------------------------------------


long	SpeciesList::NewSpecies( SpeciesList **result, short newlevel, long maxevals ) {

	SpeciesList	*head=NULL, *tmp;
	long		oldlevel = level, evals = 0,
			maxev1 = 0, maxev2 = 0;

	if( level == newlevel-1 ) // the old way if lowest level
		return _NewSpecies( result, newlevel, maxevals );

	if( maxevals < 6 ) maxev1 = maxevals; // maxev < 3 makes no sense
	else maxev1 = maxev2 = maxevals/2;
	
	// --- creating species at lowest possible level
	level = newlevel-1;
	evals += _NewSpecies( result, newlevel, maxev1 );
	if( Fail() ) return evals;

	// --- creating species at original level
	level = oldlevel;
	evals += _NewSpecies( &head, newlevel, maxev2 );
	if( Fail() ) return evals;

	// --- concatenating
	if( *result == NULL ) *result = head;
	else
	{
		for( tmp=*result; tmp->next!=NULL; tmp = tmp->next );
		tmp->next = head;
		if( head != NULL ) head->prev = tmp;
	};

	return evals;
};


// -----------------------------------------------------------------------


long	SpeciesList::_NewSpecies( SpeciesList **result, short newlevel, long maxevals ) {
// --- length is counted so that approximately
// ---     length + (length*(length-1))/2 = maxevals
// --- where (length*(length-1))/2 is the number of pairs in length elements

	SearchSpElement	**base, *oldcenter, *tmp;
	long		length, evals = 0, i, j;
	char		*toresult;
	SpeciesList	reshead,	// temporary head for result list
			*res;

	*result = NULL;
	
	FailFlag = 1==0;
	if( maxevals <= 0 ) return 0;

	FailFlag = 1==1;
	
	// --- initializing base list with random elements from attr. area
	length = (long) ( ( sqrt( 1.0 + 8*maxevals ) - 1.0 ) / 2.0 );
	oldcenter = center;
	
	base = new SearchSpElement*[ length ];
	toresult = new char[ length ];
	i = 0;

	if( base != NULL )
		for(; i<length; ++i )
		{	base[i] = INI.Prototype()->RandNew( level, oldcenter );
			if( base[i] == NULL ) break;
			++evals;
			if( base[i]->CurrValue() > center->CurrValue() ){
				center->UpdateFrom( base[i] );
			}
		};
	if( i != length )
	{
		message((char*)"No memory (NewSpecies,1)",MSG_ERROR);
		goto CLEAN;
	};
	// --- creating new species -------------------------------------
	for( i=0; i<length; ++i ) toresult[i] = 1==0;
	for( i=0; i<length-1; ++i )
	for( j=i+1; j<length; ++j )
	{
		if( toresult[i] && toresult[j] ) continue;
		tmp = INI.Prototype()->BetweenNew( base[i], base[j] );
		if( tmp == NULL )
		{
			message((char*)"No memory (NewSpecies,2)",MSG_ERROR);
			goto CLEAN;
		};
		++evals;
		if( tmp->CurrValue() > center->CurrValue() ){
			center->UpdateFrom( tmp );
		}
		else if( tmp->CurrValue() < base[i]->CurrValue() &&
		     	 tmp->CurrValue() < base[j]->CurrValue() )
			// --- new species -------------------------------
			{
				toresult[i] = toresult[j] = 1==1;
				
			}
		delete tmp;
	};

	// --- collecting result ---------------------------------------
	res = &reshead;
	for( i=0; i<length; ++i )
	{
		if( toresult[i] )
		{
			res->next = new SpeciesList( base[i], newlevel );
			if( res->next == NULL )
			{
				message((char*)"No memory (NewSpecies,3)",MSG_ERROR);
				for(; i<length; ++i ) toresult[i] = 1==0;
				goto CLEAN;
			};
			res->next->prev = res;
			res = res->next;
		};
	};
	*result = reshead.next;
	if( *result != NULL ) (*result)->prev = NULL;
	reshead.next = NULL; // to prevent result from destructing

	FailFlag = 1==0;

CLEAN:	if( base != NULL )
	{
		for( i=0; i<length && base[i]!=NULL; ++i )
			if( !toresult[i] ) delete base[i];
		delete base;
	};
	if( toresult != NULL ) delete toresult;

	return evals;
};


// -----------------------------------------------------------------------

long	SpeciesList::NewSpeciesOptiParamCereb( SpeciesList **result, short newlevel, long maxevals ) {

	SpeciesList	*head=NULL, *tmp;
	long		oldlevel = level, evals = 0,
			maxev1 = 0, maxev2 = 0;

	if( level == newlevel-1 ) // the old way if lowest level
		return _NewSpeciesOptiParamCereb( result, newlevel, maxevals );

	if( maxevals < 6 ) maxev1 = maxevals; // maxev < 3 makes no sense
	else maxev1 = maxev2 = maxevals/2;
	
	// --- creating species at lowest possible level
	level = newlevel-1;
	evals += _NewSpeciesOptiParamCereb( result, newlevel, maxev1 );
	if( Fail() ) return evals;

	// --- creating species at original level
	level = oldlevel;
	evals += _NewSpeciesOptiParamCereb( &head, newlevel, maxev2 );
	if( Fail() ) return evals;

	// --- concatenating
	if( *result == NULL ) *result = head;
	else
	{
		for( tmp=*result; tmp->next!=NULL; tmp = tmp->next );
		tmp->next = head;
		if( head != NULL ) head->prev = tmp;
	};

	return evals;
};
// -----------------------------------------------------------------------


long	SpeciesList::_NewSpeciesOptiParamCereb( SpeciesList **result, short newlevel, long maxevals ) {
// --- length is counted so that approximately
// ---     length + (length*(length-1))/2 = maxevals
// --- where (length*(length-1))/2 is the number of pairs in length elements

	SearchSpElement	**base, *oldcenter, *tmp;
	long		length, evals = 0, i, j;
	char		*toresult;
	SpeciesList	reshead,	// temporary head for result list
			*res;

	*result = NULL;
	
	FailFlag = 1==0;
	if( maxevals <= 0 ) return 0;

	FailFlag = 1==1;
	
	// --- initializing base list with random elements from attr. area
	//length = (long) ( ( sqrt( 1.0 + 8*maxevals ) - 1.0 ) / 2.0 );
	length = (long)maxevals ;
	oldcenter = center;
	base = new SearchSpElement*[ length ];
	toresult = new char[ length ];
	i = 0;
	if( base != NULL )
		for(; i<length; ++i )
		{
			base[i] = INI.Prototype()->RandNewParal( level, oldcenter );
			if( base[i] == NULL ) break;
			++evals;
			//if( base[i]->CurrValue() > center->CurrValue() )
			//	center->UpdateFrom( base[i] );
		};
	if( i != length )
	{
		message((char*)"No memory (NewSpecies,1)",MSG_ERROR);
		goto CLEAN;
	};

	// --- all the new species are inserted in the list -------------------------------------
	for( i=0; i<length; ++i ) toresult[i] = 1==1;
	/*for( i=0; i<length-1; ++i )
	for( j=i+1; j<length; ++j )
	{
		if( toresult[i] && toresult[j] ) continue;
		tmp = INI.Prototype()->BetweenNewParal( base[i], base[j] );
		if( tmp == NULL )
		{
			message((char*)"No memory (NewSpecies,2)",MSG_ERROR);
			goto CLEAN;
		};
		++evals;
		//if( tmp->CurrValue() > center->CurrValue() )
		//	center->UpdateFrom( tmp );
		//else 
		
		if( tmp->CurrValue() > base[i]->CurrValue() &&
		    tmp->CurrValue() > base[j]->CurrValue() )
		{
			base[i]->UpdateFrom( tmp );
			toresult[i] = 1==1;
		
		}
		else{	
			if( tmp->CurrValue() < base[i]->CurrValue() &&
			 tmp->CurrValue() < base[j]->CurrValue() )
				// --- new species -------------------------------
				toresult[i] = toresult[j] = 1==1;
			
		}
		
			
		delete tmp;
	};*/

	// --- collecting result ---------------------------------------
	res = &reshead;
	for( i=0; i<length; ++i )
	{
		if( toresult[i] )
		{
			res->next = new SpeciesList( base[i], newlevel); //Inicialmente seguidor asociado=lider
			if( res->next == NULL )
			{
				message((char*)"No memory (NewSpecies,3)",MSG_ERROR);
				for(; i<length; ++i ) toresult[i] = 1==0;
				goto CLEAN;
			};
			res->next->prev = res;
			res = res->next;
		};
	};
	*result = reshead.next;
	if( *result != NULL ) (*result)->prev = NULL;
	reshead.next = NULL; // to prevent result from destructing

	FailFlag = 1==0;

CLEAN:	if( base != NULL )
	{
		for( i=0; i<length && base[i]!=NULL; ++i )
			if( !toresult[i] ) delete base[i];
		delete base;
	};
	if( toresult != NULL ) delete toresult;

	return evals;
};

// -----------------------------------------------------------------------


ofstream&	SpeciesList::Save( ofstream & myfile ) {

	this->center->Save(myfile);

	myfile << "\t" << level;
};

SpeciesList*	SpeciesList::LoadFromFile(ifstream & file){

	long l;

	SearchSpElement * newCenter = NDimRealElement::LoadFromFile(file);

	file >> l;

	if (!file){
		message((char*)"Error reading SpeciesList.",MSG_ERROR);
	}

	return new SpeciesList(newCenter,l);
}


// -----------------------------------------------------------------------


void	SpeciesList::Save( FILE* stream ) {

	//FailFlag = fprintf( stream, "# level: %d\n", level ) == EOF;
	if( !FailFlag )
	{
		center->Save( stream );
		if( !center->Fail() ) return; // FailFlag is properly set!
	}

	FailFlag = 1==1;
	message( (char*)"Error writing species to stream.",MSG_ERROR);
	return;
};

// -----------------------------------------------------------------------


void	SpeciesList::Save2( FILE* stream ) {

	//FailFlag = fprintf( stream, "# level: %d\n", level ) == EOF;
	if( !FailFlag )
	{
		center->Save2( stream );
		if( !center->Fail() ) return; // FailFlag is properly set!
	}

	FailFlag = 1==1;
	message((char*) "Error writing species to stream.",MSG_ERROR);
	return;
};


// -----------------------------------------------------------------------


void	SpeciesList::Absorb( SpeciesList *sl ) {

	if( sl->center->CurrValue() > center->CurrValue() )
		center->UpdateFrom( sl->center );
};
