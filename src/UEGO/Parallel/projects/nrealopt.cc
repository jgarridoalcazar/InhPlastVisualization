#include "uego.h"

long	NDimRealElement::Optimize( short radind, long maxeval ) {
	const long	Scnt = 5, 
			Fcnt = 3;
	const double	ct = .5, 
			ex = 2.0,
			//-------------------------------------------
			//----NORMALIZAMOS RADIO!!!!!!!!!!!!!!!!!!!!!
			//-------------------------------------------
			rad=INI.R(radind)/INI.R(0)*sqrt((dim/3)*2), //normalized rad.
			sigmaub = rad,
			sigmalb = ( rad/1000.0 > 1e-5 ? 1e-5 : rad/1000.0 );

	NDimRealElement	*tmp;
	long		j, i, scnt, fcnt;
	double		sigma, *b, *epsi;

	FailFlag = 1==1;

	// --- memory allocations ----------------------------
	tmp = new NDimRealElement( dim );
	b = new double [ 2*dim ];
	if( tmp == NULL || tmp->Fail() || b == NULL )
	{
		message((char*)"No memory in NDimRealElement::Optimize()",MSG_ERROR);
		if( tmp != NULL ) delete tmp;
		return 0;
	};
	epsi = b + dim;

	scnt = 0;
	fcnt = 0;
	sigma = sigmaub;

	for( i=0; i<dim; i++ ) {
		b[i] = 0.0;
		epsi[i] = 0.0;
	}
	j = 0;
	while( (j < maxeval) && (fcnt <= 20)    ) //  && (sigma < sigmalb )
	{
		//printf("j= %ld fcnt= %ld \n", j, fcnt);


		if( scnt > Scnt ) sigma = ex * sigma;
		if( fcnt > Fcnt ) sigma = ct * sigma;
		if( sigma < sigmalb || fcnt > 20 ) // restart search
		{
			//scnt = 0;
			//fcnt = 0;
			sigma = sigmaub;
			for( i=0; i< dim; i++ ) b[i] = 0.0;
		}
		if( sigma > sigmaub ) sigma = sigmaub;
	//	printf("\t rad=%f\t sigma=%f, j<maxeval=%d\n" ,rad ,sigma,j );

		for( i=dim/3; i<dim; i++ )
		{
			epsi[i] = Gauss( b[i], sigma );
			if( epsi[i] > rad ) epsi[i] = rad;
			else if( epsi[i] < -rad ) epsi[i] = -rad;
		};

		// --- try tmp = center + epsi --------------------------
		tmp->UpdateFrom( this );
		tmp->Add( epsi );
		tmp->UpdateValue();
		j++;
		if( tmp->CurrValue() >= CurrValue() )
		{
 			UpdateFrom( tmp );
			for( i=dim/3; i < dim; ++i ) 
				b[i] = 0.4 * epsi[i] + 0.2 * b[i];
			++scnt;
			fcnt=0;
			if ( j >= maxeval ) break;
		}
		else
		{
			// --- try tmp = center - epsi --------------------
			tmp->UpdateFrom( this );
			tmp->Add( epsi, -1 );
			tmp->UpdateValue();
			j++;
			if( tmp->CurrValue() >= CurrValue() )
			{
 				UpdateFrom( tmp );
				for( i=dim/3; i < dim; ++i )
					b[i] = b[i] - 0.4*epsi[i];
				++scnt;
				fcnt=0;
			}
			else
				{
			// --- center is still better than tmp -------------
				for( i=dim/3; i < dim; ++i ) b[i] = 0.5 * b[i];
				fcnt = fcnt+1;
				scnt = 0;
				}
		}
	};

	delete tmp;
	delete b;
	FailFlag = 1==0;

	return j;
};

SearchSpElement *	NDimRealElement::InitOptimizeParal( short radind) {
	//-------------------------------------------
	//----NORMALIZAMOS RADIO!!!!!!!!!!!!!!!!!!!!!
	//-------------------------------------------
	rad=INI.R(radind); //normalized rad.
	sigmaub = rad;
	sigmalb = ( rad/1000.0 > 1e-5 ? 1e-5 : rad/1000.0 );

	FailFlag = 1==1;

	NDimRealElement * tmp = new NDimRealElement(dim);

	if (b==NULL){
		b = new double [ 2*dim ];
	}

	if( tmp == NULL || tmp->Fail() || b == NULL ){
		message((char*)"No memory in NDimRealElement::Optimize()",MSG_ERROR);
		if( tmp != NULL ) delete tmp;
		return NULL;
	};

	epsi = b + dim;
	scnt = 0;
	fcnt = 0;
	sigma = sigmaub;

	for(unsigned int i=0; i<dim; i++ ) {
		b[i] = 0.0;
		epsi[i] = 0.0;
	}

	for(unsigned int i=0; i<dim; i++ ){
		epsi[i] = Gauss( b[i], sigma );
		if( epsi[i] > rad ) epsi[i] = rad;
		else if( epsi[i] < -rad ) epsi[i] = -rad;
	};

	// --- try tmp = center + epsi --------------------------
	tmp->UpdateFrom( this );
	tmp->Add( epsi );

	this->sign = 1;

	FailFlag = 1==0;

	return tmp;
};

SearchSpElement *	NDimRealElement::ResumeOptimize( SearchSpElement *  tmpaux) {

	NDimRealElement * tmp = (NDimRealElement *) tmpaux;

	if (this->sign==1){
		if( tmp->CurrValue() >= CurrValue() ){
			this->UpdateFrom( tmp );
			for( unsigned int i=0; i < dim; ++i )
				b[i] = 0.4 * epsi[i] + 0.2 * b[i];
			++scnt;
			fcnt=0;
		} else {
			// --- try tmp = center - epsi --------------------
			this->sign = -1;
			tmp->UpdateFrom( this );
			tmp->Add( epsi, -1 );
			return tmp;
		}
	} else {
		if( tmp->CurrValue() >= CurrValue() ){
			UpdateFrom( tmp );
			for( unsigned int i=0; i < dim; ++i )
				b[i] = b[i] - 0.4*epsi[i];
			++scnt;
			fcnt=0;
		} else {
			// --- center is still better than tmp -------------
			for( unsigned int i=0; i < dim; ++i )
				b[i] = 0.5 * b[i];
			fcnt = fcnt+1;
			scnt = 0;
		}

	}

	if (fcnt>20) return NULL;

	if( scnt > Scnt ) sigma = ex * sigma;
	if( fcnt > Fcnt ) sigma = ct * sigma;
	if( sigma < sigmalb || fcnt > 20) { // restart search
		//scnt = 0;
		//fcnt = 0;
		sigma = sigmaub;
		for(unsigned int  i=0; i< dim; i++ )
			b[i] = 0.0;
	}
	if( sigma > sigmaub ) sigma = sigmaub;
	for(unsigned int  i=0; i<dim; i++ ){
		epsi[i] = Gauss( b[i], sigma );
		if( epsi[i] > rad ) epsi[i] = rad;
		else if( epsi[i] < -rad ) epsi[i] = -rad;
	};

	// --- try tmp = center + epsi --------------------------
	tmp->UpdateFrom( this );
	tmp->Add( epsi );
	return tmp;

};






