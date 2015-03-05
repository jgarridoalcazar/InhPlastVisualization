#include "uego.h"

long	NDimRealElement::Optimize( short radind, long maxeval ) {
	//printf("maxeval %ld \n", maxeval);
	//getchar();


	/*printf("OPTIMIZE ENTRADA\n");
	printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  Obj= %lf %ld\n",x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],x[20],x[21],x[22],x[23],x[24],x[25],x[26],x[27],x[28],x[29],CurrValue());
*/

/*x[0]=1.000000; x[1]=1.000000; x[2]=1.000000; x[3]=1.000000; x[4]=1.000000; x[5]=1.000000; x[6]=1.000000; x[7]=1.000000; x[8]=1.000000; x[9]=1.000000;  
x[10]=0.999959; x[11]=0.501757; x[12]=0.657321; x[13]=0.000000;  x[14]=0.467670; x[15]=0.853677; x[16]=0.591783; x[17]=0.637611; x[18]=0.000000; x[19]=0.514010;  x[20]=0.584623; x[21]=0.403368; x[22]=0.214822; x[23]=0.701864; x[24]=0.143337; x[25]=0.738816; x[26]=0.806279; x[27]=1.000000; x[28]=0.999655; x[29]=0.091454;


this->UpdateValue();


	printf("OPTIMIZE convertido\n");
	printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  Obj= %lf %ld\n",x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],x[20],x[21],x[22],x[23],x[24],x[25],x[26],x[27],x[28],x[29],CurrValue());
	
*/

//maxeval=160;
//radind = 4;

//getchar();	

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



	/*printf("OPTIMIZE SALIDA\n");
	printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  Obj= %lf %ld\n",x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],x[20],x[21],x[22],x[23],x[24],x[25],x[26],x[27],x[28],x[29],CurrValue());*/
//getchar();

	//printf("rad=%f, level=%d, INI.R(0)=%f, maxeval=%d, j=%d\n",rad,radind,INI.R(0),maxeval,j);
	return j;
	//return maxeval;
};



