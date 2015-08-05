#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include "uego.h"
#include "communication_manager.h"

// -----------------------------------------------------------------------


void	Master::ReInitParal( Ini* ini, char* trace, char *file) {

	SearchSpElement	*root;

	char		msg[250];	// working


	FailFlag = 1==1;

	// Start the communication manager thread
	to_simulate_queue.clear();
	simulated_queue.clear();
	end_simulation = false;
	if (ComManager!=0){
		ComManager->detach();
		delete ComManager;
	}
	ComManager = new CommunicationManager(to_simulate_queue, simulated_queue, end_simulation);
	ComManager->start();

	_ini = ini;		// only storing pointer!
	tracename = trace;	// only storing pointer!


	if (file==NULL){
		// --- initializing species list with root:
		// --- 1 random element from space
		head = new SpeciesList;

		if( head != NULL ) root = _ini->Prototype()->RandNewParal();

		//sprintf( msg, "MASTER :: Adding specie to the simulation queue.");
		//message( msg, MSG_INFORMATION );

		// Add the new species to the simulation list and wait
		to_simulate_queue.add(root);
		// Wait until all the simulations are finished
		root = simulated_queue.remove();

		//sprintf( msg, "MASTER :: Received specie after simulating.");
		//message( msg, MSG_INFORMATION );

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
		level = 0;

		this->generated_species = true;

		level = 0;
	} else {
		std::ifstream myFile;
		myFile.open(file);

		if (!myFile){
			message((char*)"Error on opening state file",MSG_ERROR);
			return;
		}

		myFile >> FunctionEvals >> length >> level >> generated_species;
		if (!myFile){
			message((char*)"Error on loading state file",MSG_ERROR);
			return;
		}

		// --- initializing species list from the file:
		SpeciesList	*tmp;
		head = new SpeciesList;
		if( head != NULL ) {
			tmp = head;
			for (int i=0; i<length; ++i, tmp=tmp->next){
				tmp->next = SpeciesList::LoadFromFile(myFile);
				tmp->next->prev = tmp;
				tmp->next->next = NULL;
			}
		}

		myFile.close();
		if (!myFile){
			message((char*)"Error on closing state file",MSG_ERROR);
			return;
		}
	}


	FailFlag = 1==0;
};

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//	1.CREATION
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// -----------------------------------------------------------------------


void	Master::_NewSpeciesOptiParamCereb( long evals, char *file ) {

	SpeciesList	newspechead,	// First element is newspechead.next!
			*newlst,
			*newend,	// shows end of new spec list
			*tmp;		// position in old list
	long 		newevals;

	FailFlag = 1==1;

	double * x = new double [this->_ini->Dimension()];
	char msg[25000];

	newend = &newspechead;
	for( tmp = head; tmp->next != NULL; tmp = tmp->next )
	{
		// --- to allow multiple calls w.o. level increased
		if( tmp->next->level == level ) continue;


		// --- list of new species in newlst
		newevals = tmp->next->NewSpeciesOptiParamCereb( &newlst, level, evals );
		//FunctionEvals +=newevals;
		if( tmp->next->Fail() ) return;
		
		//---Check the best
		//---El centro se compara con todos los puntos, con lo que siempre tiene el mejor.
		/*if( tmp->next->center->CurrValue() > BestObj ) {
			BestObj = tmp->next->center->CurrValue();
			BestX = tmp->next;
			SaveBest(file);
		}*/

		sprintf(msg,"New species. Individual with ");
		tmp->next->center->GetX(x);
		for (unsigned int i=0; i<this->_ini->Dimension(); ++i){
			double value = x[i];
			sprintf(msg,"%s%e\t",msg,value);
		}

		sprintf(msg,"%s%f\t%d\n",msg,tmp->next->center->CurrValue(), tmp->next->level);


		// --- insert new list to newspechead, increase length
		newend->next = newlst;
		if( newlst != NULL ) newlst->prev = newend;
		while( newend->next != NULL )
		{

			sprintf(msg,"%s generates individual with ",msg);

			newend->next->center->GetX(x);
			for (unsigned int i=0; i<this->_ini->Dimension(); ++i){
				double value = x[i];
				sprintf(msg,"%s%e\t",msg,value);
			}

			sprintf(msg,"%s. Distance: %e\n",msg,tmp->next->center->Distance(tmp->next->center,newend->next->center));

			++length;
			FunctionEvals ++;
			newend = newend->next;
		};

		printf("%s",msg);

	};

	delete [] x;

	tmp->next = newspechead.next;
	if( newspechead.next != NULL ) newspechead.next->prev = tmp;
	newspechead.next = NULL; // to prevent destructing;

	FailFlag = 1==0;
};

//--------------------------------------------------------------------------

void	Master::NewSpeciesOptiParamCereb(char *file) {

	const long	oldlength = length;
	SpeciesList *tmp;
	
	std::cout << "Evaluations: " << _ini->NewSpecEvals(level) << " Length: " << length << std::endl;

	_NewSpeciesOptiParamCereb( _ini->NewSpecEvals(level) / length ,file);
	
	if( Fail() ) return;
};



//--------------------------------------------------------------------------

void	Master::NewSpeciesParal(char *file) {
	SpeciesList     *tmp, *end;

	long	totalSpSend=0, totalSpRecv=0;
	char		msg[250];	// working
	
	int myid = INI.getmyidIni();

	sprintf( msg, "MASTER :: Creating species. Level %d length %ld",level, length );
	message( msg, MSG_INFORMATION );


	NewSpeciesOptiParamCereb(file);

	end = head->next;

	//-- reparto de especies
	totalSpSend = 0;
	for(tmp = head->next; tmp!=NULL; tmp = tmp->next)
	{
		//-- Solo las especies que se han creado en este nivel
		//-- son las que tengo que repartir --> OJO!!
		//-- Avanzamos en la lista con end, para saber la utima especie
		if( tmp->level != level ){ end=end->next; continue;}

		//sprintf( msg, "MASTER :: Adding specie to the simulation queue.");
		//message( msg, MSG_INFORMATION );

		//-- envio
		this->to_simulate_queue.add(tmp->center);
		totalSpSend++;
	};
			
	//-- ajustamos el valor objetivo a las que me quedan
	//-- y las engancho a la lista con end

	//printf(" Master :: enviadas especies correctamente \n");
	totalSpRecv = 0;
	//-- recibo la lista de los slaves
	while (totalSpSend != totalSpRecv)//(numslaves > 0)
	{
		this->simulated_queue.remove();

		//sprintf( msg, "MASTER :: Received specie after simulating.");
		//message( msg, MSG_INFORMATION );

		totalSpRecv++;
	};

//	printf("\n MASTER: Despues de recibir todo Longitud lista= %ld::\n",length);
//	for( tmp = head; tmp->next != NULL; tmp = tmp->next ){
//		double x[2];
//		tmp->next->center->GetX(x);
//		printf("%lf %lf  Obj= %lf %d\n",x[0],x[1],tmp->next->center->CurrValue(),tmp->next->level);
//	}
			
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//	3.OPTIMIZATION
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
void	Master::OptimizeParal(char *file){

	SpeciesList     *tmp;
	
	long	max_spec_num = (long)M( level );
	long	budgetSpecies = INI.Evals(level) / length ; //---Presupuesto por especie en el nivel actual

	int myid = INI.getmyidIni();
	char		msg[25000];	// working


	sprintf( msg, "MASTER :: Optimization. Level %d length %ld. Budget species %d", level, length, budgetSpecies );
	message( msg, MSG_INFORMATION );

	//-- almacenamos la longitud actual, para la creacion sig
	_length = length;


	std::map<SearchSpElement *,SpeciesList *> OptimizationMap;
	std::map<SpeciesList *,unsigned int> RemainingEvalsMap;

	unsigned int TotalEvaluations= 0;
	SearchSpElement * newSp;
	SearchSpElement * tmpSp;

	double * x = new double [this->_ini->Dimension()];

	// Generate all the new elements to be simulated and add them to the simulation queue
	if (budgetSpecies > 0){
		for(tmp = head->next; tmp != NULL; tmp = tmp->next){
			newSp = tmp->InitializeOptimize();

			if (newSp!=NULL){
				TotalEvaluations += budgetSpecies;
				RemainingEvalsMap[tmp] = budgetSpecies;
				OptimizationMap[newSp] = tmp;

				this->to_simulate_queue.add(newSp);
			} else {
				RemainingEvalsMap[tmp] = 0;
			}
		}
	}

	// Wait until simulations are finished and generate next simulations
	while (TotalEvaluations>0){
		// Wait until a simulation is finished
		SearchSpElement * tmpSp = this->simulated_queue.remove();
		//sprintf( msg, "MASTER :: Received specie after simulating.");
		//message( msg, MSG_INFORMATION );
		TotalEvaluations--;
		tmp = OptimizationMap[tmpSp];
		RemainingEvalsMap[tmp]--;

		sprintf(msg,"Optimizing. Individual with ");
		tmp->center->GetX(x);
		for (unsigned int i=0; i<this->_ini->Dimension(); ++i){
			double value = x[i];
			sprintf(msg,"%s%e\t",msg,value);
		}

		sprintf(msg,"%s%f\t%d\n",msg,tmp->center->CurrValue(), tmp->level);

		sprintf(msg,"%s generates individual with ",msg);

		tmpSp->GetX(x);
		for (unsigned int i=0; i<this->_ini->Dimension(); ++i){
			double value = x[i];
			sprintf(msg,"%s%e\t",msg,value);
		}

		sprintf(msg,"%s%f. Distance: %e\n",msg,tmpSp->CurrValue(),tmp->center->Distance(tmp->center,tmpSp));

		printf("%s",msg);

		// Generate a new individual
		newSp = tmp->ResumeOptimize(tmpSp);
		if (RemainingEvalsMap[tmp]>0){ // Check if there are simulations left and a new element has been generated
			if (newSp!=NULL){
				OptimizationMap[newSp] = tmp;

				//sprintf( msg, "MASTER :: Adding specie to the simulation queue.");
				//message( msg, MSG_INFORMATION );

				this->to_simulate_queue.add(newSp);

			} else {
				TotalEvaluations -= RemainingEvalsMap[tmp];
			}
		}
	}

	delete [] x;

	return;

};

void	Master::FinalizeParal(){

	this->end_simulation = true;

	char msg [150];
	sprintf( msg, "MASTER :: Sending ending signals to the workers");
	message( msg, MSG_INFORMATION );

	this->ComManager->join();

	sprintf( msg, "MASTER :: Ending communication manager");
	message( msg, MSG_INFORMATION );

	this->ComManager->detach();

	return;

};


//----------------------------------------------------------------------------------

void Master::packingSingleSpecies(int tag, SpeciesList * list,int tid)
{

	
	double *tmpx,tmpvalue;
	int 	i;
	short   tmplevel;
	int     position,
	        Nbytes;
	long	dim;
	char   *buf;
	
	dim = INI.Dimension();
	
	tmpx=new double[dim];
	list->center->GetX(tmpx);
	tmpvalue=list->center->CurrValue();
	tmplevel=list->level;
	
//	Nbytes=sizeof(int)+dim*sizeof(double)+sizeof(double)+sizeof(short);	
	Nbytes=dim*sizeof(double)+sizeof(double)+sizeof(short);
	buf=new char[Nbytes];
	position=0;
	
	//MPI_Pack(&dim,1,MPI_INT,buf,Nbytes,&position,MPI_COMM_WORLD);		/* Pack dim*/
	MPI_Pack(tmpx,dim,MPI_DOUBLE,buf,Nbytes,&position,MPI_COMM_WORLD);	/* Pack tmpx[dim]*/
	MPI_Pack(&tmpvalue,1,MPI_DOUBLE,buf,Nbytes,&position,MPI_COMM_WORLD);	/* Pack tmpvalue*/			
	MPI_Pack(&tmplevel,1,MPI_SHORT,buf,Nbytes,&position,MPI_COMM_WORLD);	/* Pack tmplevel*/
	MPI_Send(buf,position,MPI_PACKED,tid,tag, MPI_COMM_WORLD);		/* Send packet */
	
}

// -----------------------------------------------------------------------

void Master::sendList(){


	double		*tmpx;
  	short 		tmplevel;
	int 		ii,position,*sizelist,*Nbytes;
	char		*sendBuf;
	SpeciesList	*tmp;
	MPI_Status	status;
	int		m,n,dim,auxSizeList,remainder;
	
	
	
	m = INI.getmyidIni();
	n = INI.getnumprocIni();  
	dim = INI.Dimension(); 
	
	tmpx=new double[dim];
	sizelist = new int[n];
	Nbytes = new int[n];

	//---Determine the total size of the list
	auxSizeList = 0;
        for( tmp = head; tmp->next!= NULL; tmp = tmp->next )    auxSizeList++;
	
	printf(" Total sizelist %i\n", auxSizeList);
	remainder = auxSizeList%n; 		
	auxSizeList = auxSizeList/n;
        	

        //---Update values of sizelist for each processor
        for (int ii = 0; ii < n; ii++) sizelist[ii] = auxSizeList;
        for (int ii=0 ; ii < remainder  ; ii++) sizelist[ii] ++;
	
	
	printf(" Partial sizelist %i\n", auxSizeList);
	printf(" Remainder %i\n", remainder);
	for (int ii = 0; ii < n; ii++){
		printf(" Tamaño enviado al proc %i es de %i\n", ii,sizelist[ii] );
		Nbytes[ii] = sizelist[ii]*(dim*sizeof(double)+sizeof(double)); //---Compute de number of bytes for each processor
		MPI_Send(&sizelist[ii],1,MPI_INT,ii,56,MPI_COMM_WORLD);	   //---Send the length of the sublist to each processor
	}
	
	


	/*tmp=head;
        	for (long ii = 0; ii < n; ii++){
            parameterInput[ii].tmphead = tmp->next;
            for(long aa = 0; aa < parameterInput[ii].sizelist; aa++) tmp=tmp->next;



	Nbytes=size*(dim*sizeof(double)+sizeof(double));	
	tmp = head;
	
	for(tmp=head,ii=tid*size;tmp->next!=NULL;tmp=tmp->next,sizelist++); 
	
	for(tmp=head,sizelist=0;tmp->next!=NULL;tmp=tmp->next,sizelist++); 
		//printf("sizelist= %ld\n", sizelist);
		Nbytes=sizelist*(dim*sizeof(double)+sizeof(short));	
		Nbytes2=2*sizeof(int);	
	
		sendBuf=(char*)malloc(Nbytes*sizeof(char));
		sendBuf2=(char*)malloc(Nbytes2*sizeof(char));	
	
		position=0;
		MPI_Pack(&m,1,MPI_INT,sendBuf2,Nbytes2,&position,MPI_COMM_WORLD);		
		MPI_Pack(&sizelist,1,MPI_INT,sendBuf2,Nbytes2,&position,MPI_COMM_WORLD);		
		MPI_Send(sendBuf2,position,MPI_PACKED,0,27,MPI_COMM_WORLD);
	
		
		position=0;	
		for(tmp=head->next;tmp!=NULL;tmp=tmp->next){
			tmp->center->GetX(tmpx);
			//for(long aa=0;aa<dim;aa++) tmpx[aa] = tmpx[aa]*(INI.Upb(aa)-INI.Lowb(aa)) + INI.Lowb(aa); // --- Normalizo antes de enviar
			tmplevel=tmp->level;
			MPI_Pack(tmpx,dim,MPI_DOUBLE,sendBuf,Nbytes,&position,MPI_COMM_WORLD);	
			MPI_Pack(&tmplevel,1,MPI_SHORT,sendBuf,Nbytes,&position,MPI_COMM_WORLD);
 		}	
		if(MPI_Send(sendBuf,position,MPI_PACKED,0,26,MPI_COMM_WORLD)!=MPI_SUCCESS)
			printf("An error in MPI_Send() sincronizacion\n");
	
		// ---Free memory
		free(tmpx);
		free(sendBuf);
		free(sendBuf2);
	
		// ---Elimino la lista que acabo de enviar.	
		if( head != NULL ){
			while( head->next != NULL ){
				tmp = head->next;
				head->next = head->next->next;
				delete tmp;
			};
		};*/
	
	
}//end sendOptim



