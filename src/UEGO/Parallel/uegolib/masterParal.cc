#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "uego.h"

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


		// --- insert new list to newspechead, increase length
		newend->next = newlst;
		if( newlst != NULL ) newlst->prev = newend;
		while( newend->next != NULL )
		{
			++length;
			FunctionEvals ++;
			newend = newend->next;
			
		};
	};
	tmp->next = newspechead.next;
	if( newspechead.next != NULL ) newspechead.next->prev = tmp;
	newspechead.next = NULL; // to prevent destructing;

	FailFlag = 1==0;
};

//--------------------------------------------------------------------------

void	Master::NewSpeciesOptiParamCereb(char *file) {

	const long	oldlength = length;
	SpeciesList *tmp;
	
	_NewSpeciesOptiParamCereb( _ini->NewSpecEvals(level) / length ,file);
	
	if( Fail() ) return;
};



//--------------------------------------------------------------------------

void	Master::NewSpeciesParal(char *file) {
	SpeciesList     *tmp;

	//------------------------------------------
	//------ Parallel variables for MPI --------
	//------------------------------------------
	int	namelen=0, myid=0, numproc=0, numslaves=0,
		tag_level=1, i=0,j=0, turn=0;
			
	char	name[25];
	
	char	 *sendBuf=NULL,*send=NULL,
		 *recvBuf=NULL,*recv=NULL;;
	
	long	totalSpSend=0, totalSpRecv=0,total_followerGo=0;
		
	MPI_Status	 status;
	
	double	timer1=0.0, timer2=0.0, timer3=0.0,
		value=0.0;//,	radio=0.0;
	double	*tmpx = NULL;
	int	tmplevel=0, tmpfevals=0, id_recv=0, nSp_master=0,
		Nbytes = 0,
		position=0, correct=0,	fin=0,		
		dim=3, cont=0, n=0, remainder=0,
		numSpecies=0, numSpRecv=0;
	
	SpeciesList	*newlst=NULL, *newend=NULL, newspechead,
			*act=NULL, *tmp2=NULL, *end=NULL;
	
	SearchSpElement	*newSp_center=NULL, *newSp_follower=NULL;
	

	SearchSpElement *opti_center, *ini_center;
	double		 ini_value;
	char		msg[250];	// working
	long	newevals,
		max_spec_num = (long)M( level );	
	long	bugetSpecies = INI.Evals(level) / max_spec_num ; //---Presupuesto por especie en el nivel actual
	
	//---Obtaining some important values
	myid = INI.getmyidIni();
	numproc = INI.getnumprocIni();  
	dim = INI.Dimension(); 
	

	//-- Allocation memory
	tmpx=new double[dim];
	Nbytes = (dim*sizeof(double)) + sizeof(short) +sizeof(double);
	sendBuf=(char*)malloc(Nbytes*sizeof(char));
	recvBuf=(char*)malloc(Nbytes*sizeof(char));
		
	
	printf(" Hola \n");


			

	if(myid == 0)
	{	
		//printf("****************************************************\n");
				
		//------------------------
		//-------- MASTER --------
		//------------------------
		sprintf( msg, "P %2d MASTER :: Creating species. Level %d length %ld",myid, level, length );
		message( msg, MSG_INFORMATION );
				
		
		NewSpeciesOptiParamCereb(file);


		printf("\n Longitud lista= %ld::\n",length);
		for( tmp = head; tmp->next != NULL; tmp = tmp->next )
		{
			double x[2];
			tmp->next->center->GetX(x);
			printf("%lf %lf  Obj= %lf %d\n",x[0],x[1],tmp->next->center->CurrValue(),tmp->next->level);
		}
		//getchar();
		
		//************************************************************
		//step=1;
		//timer[lider->level][step][0] = MPI_Wtime();
		//*****************************	
		//lider->ini().SetFlagWho(0);
		//lider->ini().SetNewCenter(tmpxF);
		//lider->NewSpeciesLider();
		//*****************************		
		//timer[lider->level][step][1] = MPI_Wtime();
		//************************************************************
				
				
				
		//-- reparto de especies
		//************************************************************
		//step=2;
		//timer[lider->level][step][0] = MPI_Wtime();				
		totalSpSend = 0;				
		newend = &newspechead;
		cont = 0;
		numslaves = numproc-1;
		end = head->next;				
		for(tmp = head->next; tmp!=NULL; tmp = tmp->next)
		{
			//-- Solo las especies que se han creado en este nivel
			//-- son las que tengo que repartir --> OJO!!
			//-- Avanzamos en la lista con end, para saber la utima especie
			if( tmp->level != level ){ end=end->next; continue;}
					
			//-- envio 
			tmp->center->GetX(tmpx);
			value = 0.0;//tmp->center->CurrValue();
					
			if(numslaves == 0)
			{
				//-- se lo queda el master para ajustar el seguidor
				//tmp->center->UpdateValue();			//Actualizo el valor de la función objetivo
				//value = tmp->center->CurrValue();
				newSp_center=INI.Prototype()->SetNew(value,tmpx); 				
				newend->next= new SpeciesList(newSp_center, level);				
				if( newend->next != NULL ) newend->next->prev = newend;
				newend=newend->next;	
				cont++;					
				//printf("\tKeep: cont %d::Master = %.10le %.10le  %.10le	level slave %d\n",
				//	cont, tmpx[0],tmpx[1],value, level);						
			}else{
				//-- lo envio al slave "numsalves"
				fin = 0;// NO fin del bucle
				if(MPI_Send(&fin,1,MPI_LONG,numslaves,10,MPI_COMM_WORLD)!=MPI_SUCCESS)
					printf("An error in MPI_Send() sincronizacion\n");
				if(MPI_Send(tmpx, dim, MPI_DOUBLE, numslaves, 11, MPI_COMM_WORLD)!=MPI_SUCCESS)
					printf("An error in MPI_Send() send x lider\n");
				//printf("\tSend %d::Lider = %.10le %.10le %.10le\n",
				//		numslaves, tmpxL[0],tmpxL[1],tmpxL[2]);

				//printf("\nSend: cont %d::Slave = %.10le %.10le  %.10le	level slave %d\n",
				//	cont, tmpx[0],tmpx[1],value, level);
			}
			totalSpSend++;
					
			/*if (numslaves ==1) {
				if(numproc==2) numslaves=0;
				else	       numslaves=numproc-1;
			}
			else{
				if( numslaves==0) numslaves=numproc-1;
				else numslaves--;
			}*/	
			if( numslaves==0) numslaves=numproc-1;
			else numslaves--;		
		};
				
		fin=1; // fin del bucle
		for(i=1; i<numproc; i++){
			if(MPI_Send(&fin, 1, MPI_LONG, i, 10, MPI_COMM_WORLD)!=MPI_SUCCESS)
				printf("An error in MPI_Send() sincronizacion\n");
		}				
		//timer[lider->level][step][1] = MPI_Wtime();
		//************************************************************
				
				
		//************************************************************
		//step=3;
		//timer[lider->level][step][0] = MPI_Wtime();
		

		//-- ajustamos el valor objetivo a las que me quedan
		//-- y las engancho a la lista con end

		printf(" Master :: enviadas especies correctamente \n"); 
		totalSpRecv = 0;
		tmp = newspechead.next;
		while(tmp!=NULL)
		{
			//ini_lider_center = tmp->getCenter();
			//ini_lider_center->GetX_norm(tmpxL);

			tmp->center->GetX(tmpx);
			tmp->center->UpdateValue();			//Actualizo el valor de la función objetivo
			tmp->level = level;
			value = tmp->center->CurrValue();				
			//printf(" Valor objetivo = %lf \n", value);			
			tmp = tmp->next;

			//Escribimos en la lista global
			end->center = INI.Prototype()->updateCenter(tmpx, value);
			end->level = level;
			end = end->next;
			totalSpRecv++;					
		}				
		//timer[lider->level][step][1] = MPI_Wtime();
		//************************************************************
				
		printf(" Master :: evaluadas especies correctamente \n");		
			
				
		//************************************************************
		//step=4;
		//timer[lider->level][step][0] = MPI_Wtime();
		//-- recibo la lista de los slaves
		fin=0;
		numslaves = numproc-1;
		while (totalSpSend != totalSpRecv)//(numslaves > 0) 
		{														
			MPI_Iprobe(MPI_ANY_SOURCE,12,MPI_COMM_WORLD,&correct,&status);
			if(correct==1)
			{	
				id_recv = status.MPI_SOURCE;	
				MPI_Recv(&numSpRecv,1, MPI_INT, id_recv, 12, MPI_COMM_WORLD,&status); 						
				printf("\n#### Slave %d	numSpRecv %d\n", status.MPI_SOURCE, numSpRecv);
				numslaves--;
						
				Nbytes = numSpRecv * (dim*sizeof(double) + sizeof(double) );				
				recvBuf=(char*)malloc(Nbytes*sizeof(char));
				MPI_Recv(recvBuf, Nbytes, MPI_PACKED, id_recv, 13, MPI_COMM_WORLD,&status); 									
				position=0;						
				for(j=0; j<numSpRecv && end!=NULL; j++, end = end->next)
				{						
					MPI_Unpack(recvBuf,Nbytes,&position,tmpx,dim,MPI_DOUBLE,MPI_COMM_WORLD); 
					MPI_Unpack(recvBuf,Nbytes,&position,&value,1,MPI_DOUBLE,MPI_COMM_WORLD); 
			
					//printf("@@@@@RecvM::Master = %.10le %.10le %.10le  \t",
					//	tmpx[0],tmpx[1],value);
					end->center = INI.Prototype()->updateCenter(tmpx, value);
					end->level = level;
					printf("Valor objetivo end %lf \n", end->center->CurrValue());
					//getchar();
		
					totalSpRecv++;
																
				}
				printf(" Master :: recibidas especies correctamente \n");					
				//if(recvBuf!=NULL) free(recvBuf);
			}																	
		};//end_while_numslaves_fin
		//timer[lider->level][step][1] = MPI_Wtime();
		//************************************************************
				
		//--Liberar memoria
		if(totalSpRecv>0) newend=0;

		printf("\n MASTER: Despues de recibir todo Longitud lista= %ld::\n",length);
			for( tmp = head; tmp->next != NULL; tmp = tmp->next )
			{
				double x[2];
				tmp->next->center->GetX(x);
				printf("%lf %lf  Obj= %lf %d\n",x[0],x[1],tmp->next->center->CurrValue(),tmp->next->level);
			}

		//getchar();





				
	}else{
		//------------------------
		//-------- SLAVES --------
		//------------------------
				
		//************************************************************
		//step = 2;
		//timer[lider->level][step][0] = MPI_Wtime();
		//-- Recibimos especies
		fin=0;
		totalSpRecv = 0;
		newend = &newspechead;	
		do{
			MPI_Recv(&fin, 1, MPI_LONG,0, 10, MPI_COMM_WORLD, &status); 
			if(fin == 0){
				MPI_Recv(tmpx,sizeof(tmpx), MPI_DOUBLE, 0, 11, MPI_COMM_WORLD, &status);
				value = 0.0;					//Este es el valor que quiero actualizar. Por lo pronto lo pongo a cero
				newSp_center=INI.Prototype()->SetNew(value,tmpx);    
				newend->next= new SpeciesList(newSp_center, level);
				
				if( newend->next != NULL ) newend->next->prev = newend;
				newend=newend->next;
				totalSpRecv++; //-- OJO lista no esta en head->next!						
				//printf("\tRecv %2d::Lider = %.10le %.10le %.10le	level slave %d\n",
				//	myid, tmpxL[0],tmpxL[1],tmpxL[2], lider->level);
			}
		}while(fin == 0);

		printf(" Slave :: recibidas especies correctamente \n");
		//timer[lider->level][step][1] = MPI_Wtime();
		//************************************************************				
				
				
		//************************************************************
		//step=3;
		//timer[lider->level][step][0] = MPI_Wtime();
		//-- Calculamos valor objetivo y empaquetamos!!	
		if(totalSpRecv > 0)
		{							
			Nbytes = totalSpRecv * (dim*sizeof(double) + sizeof(double) );
			sendBuf=(char*)malloc(Nbytes*sizeof(char));			
			position=0;
			i=0;
			tmp = newspechead.next;
			while( (tmp != NULL) && (i < totalSpRecv) )
			{


				tmp->center->GetX(tmpx);
				tmp->center->UpdateValue();			//Actualizo el valor de la función objetivo
				value = tmp->center->CurrValue();				
			
				
				MPI_Pack(tmpx, dim, MPI_DOUBLE, sendBuf, Nbytes,&position,MPI_COMM_WORLD);
				MPI_Pack(&value, 1, MPI_DOUBLE, sendBuf, Nbytes,&position,MPI_COMM_WORLD);	
						
				printf("P %2d::Slave = %.10le %.10le %.10le level %d\tFollower = %.10le %.10le %.10le %.10le\n",myid, tmpx[0],tmpx[1],value, tmplevel);												
				tmp=tmp->next; i++;					
			};//end_while_Go
			printf(" Slave: Species empaquetas correctamente \n");
		
			if( MPI_Send(&totalSpRecv, 1, MPI_INT, 0, 12, MPI_COMM_WORLD) != MPI_SUCCESS)
				printf("An error in MPI_Send() Master :: send list to master\n");
			printf(" Slave: Total enviadas = %ld \n", totalSpRecv);
			if( MPI_Send(sendBuf, position, MPI_PACKED, 0, 13, MPI_COMM_WORLD) != MPI_SUCCESS)
				printf("An error in MPI_Send() send list to master\n");				
			printf(" Slave: Species enviadas = %ld \n", totalSpRecv);		
			//if(sendBuf!=NULL) free(sendBuf);sendBuf = NULL;	
			newend=0;
		}
		//timer[lider->level][step][1] = MPI_Wtime();
		//************************************************************	
				
				
						
	}//end_if_myid
			
			
	// --- All the proccess have to wait here.
	MPI_Barrier( MPI_COMM_WORLD );
			
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//	3.OPTIMIZATION
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
void	Master::OptimizeParal(char *file){

	SpeciesList     *tmp;

	//------------------------------------------
	//------ Parallel variables for MPI --------
	//------------------------------------------
	int	namelen=0, myid=0, numproc=0, numslaves=0,
		tag_level=1, i=0,j=0, turn=0;
			
	char	name[25];
	
	char	 *sendBufOpt=NULL,*sendBuf=NULL,
		 *recvBufOpt=NULL,*recvBuf=NULL;;
	
	long	totalSpSend=0, totalSpRecv=0,total_followerGo=0;
		
	MPI_Status	 status;
	
	double	timer1=0.0, timer2=0.0, timer3=0.0,
		value=0.0;//,	radio=0.0;
	double	*tmpx = NULL;
	int	tmplevel=0, tmpfevals=0, id_recv=0, nSp_master=0,
		NbytesOpt=0, Nbytes = 0,
		position=0, correct=0,	fin=0,		
		dim=3, cont=0, n=0, remainder=0,
		numSpecies=0, numSpRecv=0;
	
	SpeciesList	*newlst=NULL, *newend=NULL, newspechead,
			*act=NULL, *tmp2=NULL, *end=NULL;
	
	SearchSpElement	*newSp_center=NULL, *newSp_follower=NULL;
	

	SearchSpElement *opti_center, *ini_center;
	double		 ini_value;
	char		msg[250];	// working
	long	newevals,
		max_spec_num = (long)M( level );	
	long	bugetSpecies = INI.Evals(level) / max_spec_num ; //---Presupuesto por especie en el nivel actual
	
	//---Obtaining some important values
	myid = INI.getmyidIni();
	numproc = INI.getnumprocIni();  
	dim = INI.Dimension(); 
	

	//-- Allocation memory
	tmpx=new double[dim];
	NbytesOpt = (dim*sizeof(double)) + sizeof(short) +sizeof(double);
	sendBufOpt=(char*)malloc(NbytesOpt*sizeof(char));
	recvBufOpt=(char*)malloc(NbytesOpt*sizeof(char));
		
	

			
	if(myid == 0)
	{	
		//------------------------
		//-------- MASTER --------
		//------------------------
		sprintf( msg, "P %2d MASTER :: Optimization. Level %d length %ld",myid, level, length );
		message( msg, MSG_INFORMATION );
				
		//-- almacenamos la longitud actual, para la creacion sig
		_length = length;		
				
		//imprimirLista();
				
		//************************************************************
		//step=8;
		//timer[lider->level][step][0] = MPI_Wtime();
		//-- reparto de especies
		totalSpSend = 0;	
		totalSpRecv = 0; //las del master			
		newend = &newspechead;	
		cont = 0;
		numslaves = numproc-1;
		for(tmp = head->next; tmp != NULL; tmp = tmp->next)
		{
			//-- envio todo!
			tmp->center->GetX(tmpx);
			value = tmp->center->CurrValue();
			tmplevel=tmp->level;
			//Param-names is a constant! I think I do not need to send them
		
			if(numslaves == 0)
			{
				//-- se lo queda el master para ajustar el seguidor
				newSp_center=INI.Prototype()->SetNew(value,tmpx);    // Nuevo centro = tmpx1L	
				newend->next= new SpeciesList(newSp_center, tmplevel);				
				if( newend->next != NULL ) newend->next->prev = newend;
				newend=newend->next;						
				printf("\tKeep cont %d:: Solution= %.10le %.10le  %.10le  %d\n",
					cont, tmpx[0],tmpx[1],value, tmplevel);
				totalSpRecv++;	
				cont++;		
			}else{
				//-- lo envio al slave "numsalves"
				fin = 0;// NO fin del bucle
				position = 0;						
				MPI_Pack(tmpx, dim, MPI_DOUBLE, sendBufOpt, NbytesOpt,&position,MPI_COMM_WORLD);
				MPI_Pack(&value, 1, MPI_DOUBLE, sendBufOpt, NbytesOpt,&position,MPI_COMM_WORLD);	
				MPI_Pack(&tmplevel,1,MPI_SHORT, sendBufOpt, NbytesOpt,&position,MPI_COMM_WORLD);
						
				if(MPI_Send(&fin,1,MPI_LONG,numslaves,30,MPI_COMM_WORLD)!=MPI_SUCCESS)
					printf("An error in MPI_Send() sincronizacion\n");
							
				if( MPI_Send(sendBufOpt, position, MPI_PACKED, numslaves, 31, MPI_COMM_WORLD) != MPI_SUCCESS)
					printf("An error in MPI_Send() Master :: sendSpToSlaves\n");
						
				printf("\tSend %d::Solution = %.10le %.10le %.10le  %d\n",
					numslaves, tmpx[0],tmpx[1],value, tmplevel);
			}
			totalSpSend++;
					
			if( numslaves==0) numslaves=numproc-1;
			else numslaves--;		
		};//end_for_reparto
		
		//---Le indicamos a los slaves que no van a recibir mas especies		
		fin=1; // fin del bucle
		for(i=1; i<numproc; i++){
			if(MPI_Send(&fin, 1, MPI_LONG, i, 30, MPI_COMM_WORLD)!=MPI_SUCCESS)
				printf("An error in MPI_Send() sincronizacion\n");
		}
		//timer[lider->level][step][1] = MPI_Wtime();
		//************************************************************
				
				
		//-- sobre escribo la lista! ojo!
		end = head->next; 
				
		//printf(" TotalSend= %ld -- TotalReceive= %ld\n", totalSpSend, totalSpRecv);
				
		//************************************************************
		//step=9;
		//timer[lider->level][step][0] = MPI_Wtime();
		//-- Si me he quedado con alguna, las optimizo!
		if(totalSpRecv > 0)
		{

			
			tmp = newspechead.next;					
			while(tmp!=NULL)
			{
				// ---Obtengo lider y seguidor iniciales
				/*ini_center = tmp->getCenter(); 
				ini_center->GetX(tmpxL);					
				ini_value = ini_center->CurrValue();
				tmplevel = tmp->level;
				*/
				newevals=tmp->Optimize(bugetSpecies);
				FunctionEvals+=newevals;
				if( tmp->Fail() ) return;
				
				//--- Avanzamos en la lista
				tmp = tmp->next;
				end = end->next;
			}
			newend=0;
		};//end_if_totalSpRecv>0



		//timer[lider->level][step][1] = MPI_Wtime();
		//************************************************************
				
		//printf("\nLISTA FINAL DESPUES DE OPTIMIZAR MSATER:\n");
		//lider->imprimirLista();			
				
				
		//printf("@@ TotalSend= %ld -- TotalReceive= %ld\n", totalSpSend, totalSpRecv);
				
		//************************************************************
		//step=12;
		//timer[lider->level][step][0] = MPI_Wtime();
		//-- recibo la lista de los slaves
		fin=0;
		numslaves = numproc-1;
		while (totalSpSend != totalSpRecv) 
		{	
			//************************************************************
			//timer[lider->level][13][0] = MPI_Wtime();						
			MPI_Iprobe(MPI_ANY_SOURCE,40,MPI_COMM_WORLD,&correct,&status);
			if(correct==1)
			{	
				id_recv = status.MPI_SOURCE;	
				MPI_Recv(&numSpRecv, 1,MPI_INT, id_recv, 40,  MPI_COMM_WORLD,&status); 						
				printf("\n#### Slave %d	numSpRecv %d\n", status.MPI_SOURCE, numSpRecv);
				numslaves--;
						
				Nbytes = numSpRecv * (dim*sizeof(double) + sizeof(short) + sizeof(double) )+ sizeof(long)  ;
				recvBuf=(char*)malloc(Nbytes*sizeof(char));
				MPI_Recv(recvBuf, Nbytes, MPI_PACKED, id_recv, 41, MPI_COMM_WORLD, &status); 									
				position=0;					
				for(j=0; j<numSpRecv && end!=NULL; j++, end = end->next)
				{							
					MPI_Unpack(recvBuf,Nbytes,&position,tmpx,dim,MPI_DOUBLE,MPI_COMM_WORLD); 
					MPI_Unpack(recvBuf,Nbytes,&position,&value,1,MPI_DOUBLE,MPI_COMM_WORLD); 			
					MPI_Unpack(recvBuf,Nbytes,&position,&tmplevel,1,MPI_SHORT,MPI_COMM_WORLD);			
					//printf("\tRecvM::Lider = %.10le %.10le %.10le %.10le  %d\tFollower = %.10le %.10le %.10le %.10le\n",
					//	tmpxL[0],tmpxL[1],tmpxL[2],valueL, tmplevel, tmpxF[0],tmpxF[1],tmpxF[2], valueF);
					
					end->center = ini().Prototype()->SetNew(value,tmpx);
					end->level = tmplevel;
					totalSpRecv++;														
				}	
				MPI_Unpack(recvBuf,Nbytes,&position,&tmpfevals,1,MPI_LONG,MPI_COMM_WORLD); 
				FunctionEvals += tmpfevals;												
				//if(recvBuf!=NULL) free(recvBuf);						
						
			}
							
		};//end_while_numslaves_fin
		//timer[lider->level][step][1] = MPI_Wtime();
		//************************************************************
				
				
				
		//************************************************************
		//timer[lider->level][13][0] = MPI_Wtime();
		//	lider->orderListLider(lider->level);
		//timer[lider->level][13][1] = MPI_Wtime();
		//************************************************************	
		//************************************************************
		///timer[lider->level][14][0] = MPI_Wtime();
		//lider->FuseLD();
		//timer[lider->level][14][1] = MPI_Wtime();
		//************************************************************			
		//printf("LISTA GLOBAL DESPUES DE RECIBIR ESPECIES:\n");
		//lider->imprimirLista();
				
						
	}else{
		//------------------------
		//-------- SLAVES --------
		//------------------------
				
				
		//************************************************************
		//step=8;
		//timer[lider->level][step][0] = MPI_Wtime();
		//-- Recibimos especies
		fin=0;
		length = 0; //usamos como variable local
		totalSpRecv=0;
		newend = &newspechead;
		do{
			MPI_Recv(&fin, 1, MPI_LONG,0, 30, MPI_COMM_WORLD, &status); 
			if(fin == 0)
			{
				MPI_Recv(recvBufOpt,NbytesOpt,MPI_PACKED,0, 31, MPI_COMM_WORLD,&status); 	
				position=0;						
				MPI_Unpack(recvBufOpt,NbytesOpt,&position,tmpx,dim,MPI_DOUBLE,MPI_COMM_WORLD); 
				MPI_Unpack(recvBufOpt,NbytesOpt,&position,&value,1,MPI_DOUBLE,MPI_COMM_WORLD); 			
				MPI_Unpack(recvBufOpt,NbytesOpt,&position,&tmplevel,1,MPI_SHORT,MPI_COMM_WORLD);
				printf("\tRecv::Solution = %.10le %.10le  %.10le  %d\n",
					tmpx[0],tmpx[1], value, tmplevel);
				
				newSp_center=INI.Prototype()->SetNew(value,tmpx);    // Nuevo centro = tmpx1L
				
				newend->next= new SpeciesList(newSp_center, tmplevel);
				
				if( newend->next != NULL ) newend->next->prev = newend;
				newend=newend->next;
				totalSpRecv++;
			}
		}while(fin == 0);				
				



		//if(recvBuf!=NULL) free(recvBuf); recvBuf=NULL;
		//timer[lider->level][step][1] = MPI_Wtime();
		//************************************************************
				
		//-- Optimizo las especies que me he quedado
		if(totalSpRecv > 0)
		{

			printf("Optimizing... myid= %i \n", myid);
			Nbytes = totalSpRecv * (dim*sizeof(double) + sizeof(short) + sizeof(double) )+sizeof(long);
			sendBuf=(char*)malloc(Nbytes*sizeof(char));			
			position=0;
			FunctionEvals = 0;
			tmp = newspechead.next;
						
			while(tmp!=NULL)
			{
				// ---Obtengo lider y seguidor iniciales
				/*ini_center = tmp->getCenter(); 
				ini_center->GetX(tmpxL);					
				ini_value = ini_center->CurrValue();
				tmplevel = tmp->level;
				*/
				newevals=tmp->Optimize(bugetSpecies);
				FunctionEvals+=newevals;
				if( tmp->Fail() ) return;
				
				//-- empaqueto lo que tengo que enviar despues
				tmplevel=tmp->level;
				value = tmp->center->CurrValue();
				tmp->center->GetX(tmpx);			
				
				MPI_Pack(tmpx, dim, MPI_DOUBLE, sendBuf, Nbytes,&position,MPI_COMM_WORLD);
				MPI_Pack(&value, 1, MPI_DOUBLE, sendBuf, Nbytes,&position,MPI_COMM_WORLD);	
				MPI_Pack(&tmplevel,1,MPI_SHORT, sendBuf, Nbytes,&position,MPI_COMM_WORLD);
																
				printf("P %2d::Solution = %.10le %.10le %.10le level %d\n",
					myid, tmpx[0],tmpx[1], value, tmplevel);
						
				//--- Avanzamos en la lista
				tmp = tmp->next;
			};//end_while_optimize					
			MPI_Pack(&FunctionEvals, 1, MPI_LONG, sendBuf, Nbytes,&position,MPI_COMM_WORLD);
			//timer[lider->level][step][1] = MPI_Wtime();
			//************************************************************
					
					
			//************************************************************
			//step=12;
			//timer[lider->level][step][0] = MPI_Wtime();
			if( MPI_Send(&totalSpRecv, 1, MPI_INT, 0, 40, MPI_COMM_WORLD) != MPI_SUCCESS)
				printf("An error in MPI_Send() Master :: send list to master\n");
			if( MPI_Send(sendBuf, position, MPI_PACKED, 0, 41, MPI_COMM_WORLD) != MPI_SUCCESS)
				printf("An error in MPI_Send() Master :: sendSpToSlaves\n");
					
			//if(sendBuf!=NULL) free(sendBuf); sendBuf=NULL;		
			newend=0;
								
			//timer[lider->level][step][1] = MPI_Wtime();
			//************************************************************
		}//end_if_totalSpRecv>0
				
	};//end_else_myid
			
			
			
	// --- All the proccess have to wait here.
	MPI_Barrier( MPI_COMM_WORLD );




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



