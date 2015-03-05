/*----------------------------------------------------------------------------*/
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
                                                                              */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <values.h>

/* Global variables. */

#define DEBUG 1
#define ZERO 0.001
#define CELL 165                /* Dimension of the unit cell: 165 Angstroms. */

struct {
    float Dresmin, 
          Dresmax;
    int   IQmax;
    } Limits;
/*----------------------------------------------------------------------------*/
void Usage(int, char **);
void Skipline(FILE *);
void FriedelInversion(int *h, int *k, float *phs);
int ReadFile(char *infile, int **ih, int **ik, 
                           float **amp, float **phs, int **iq);
int CompReferencePHS(int Ni, int *hi, int *ki, float *ampi, float *, int *iqi,
                      int Nr, int *hr, int *kr, float *ampr, float *phsr,  
                      int *iqr, float **phsc, float **wgt, float **respot, 
                      float **iptest, float **lspec);
/*----------------------------------------------------------------------------*/
main(int argc, char **argv)
{
char ifile[100], rfile[100];
int Ni, Nr;
int Ncomp;
int *hi, *ki, *iqi;      /* H, K, IQ de la imagen a alinear. */
float *ampi, *phsi;      /* Amplitud y Fase de la imagen a alinear. */

int *hr, *kr, *iqr;      /* H, K, IQ de la imagen de referencia. */
float *ampr, *phsr;      /* Amplitud y Fase de la imagen de referencia. */

float *phsc, *wgt, *respot, *iptest, *lspec;

Limits.IQmax = 8;
Limits.Dresmin = 3.0;  /* Minimum resolution: 10 Angstroms  */
Limits.Dresmax = CELL;

/* Checks the command line. */
if(argc != 3){
    ERRHandler(0);                        /* Wrong arguments in command line. */
    Usage(argc, argv);
    exit(1);
    }

/* Name of input file: Imagen a alinear. */
strcpy(ifile, argv[1]);

/* Name of input file: Imagen de Referencia. */
strcpy(rfile, argv[2]);

/* Lectura de la imagen a alinear. */
Ni = ReadFile(ifile, &hi, &ki, &ampi, &phsi, &iqi);

/* Lectura de la imagen de referencia. */
Nr = ReadFile(rfile, &hr, &kr, &ampr, &phsr, &iqr);

/* Cálculo de PHSC: las fases para efectuar la comparación, y el resto de
   arrays auxiliares. */
Ncomp = CompReferencePHS(Ni, hi, ki, ampi, phsi, iqi,
                         Nr, hr, kr, ampr, phsr, iqr, 
                         &phsc, &wgt, &respot, &iptest, &lspec);

if(DEBUG)
    printf("\n Numero comparaciones: %5d\n", Ncomp);

/* Ending the program. */
exit(0);
}
/*----------------------------------------------------------------------------*/
void Usage(int argc, char **argv)
{
printf("\nGoal: To read two .aph files: image to align and reference image.\n");
printf("\nUsage:\n");
printf("\n\t%s ifile rfile\n", argv[0]);
printf("\n\tWhere:\n");
printf("\n\t\tifile  - Name of image to align.");
printf("\n\t\trfile  - Name of reference image.");
printf("\n\n");
}
/*----------------------------------------------------------------------------*/
int ReadFile(char *infile, int **ih, int **ik, 
                           float **amp, float **phs, int **iq)
{
char line[200];
FILE *aphf;                               /* Input file containing the spots. */
int result;                                         /* Result of the reading. */
int nlines;                                       /* No. lines in input file. */
int indx;

/* Opening the input file. */
if(!(aphf = fopen(infile, "r"))){
    ERRHandler(1);                              /* Unable to open input file. */
    exit(1);
    } 
    
    
/* How many lines? */
nlines=-1;                           /* First line is the title: Not counted. */
while(!feof(aphf)) {
 fgets(line, 200, aphf);
 nlines++;
 }
nlines--;                                            /* Last line is garbage. */

/* Allocating memory. */
*ih    = (int *)malloc(nlines * sizeof(int));
*ik    = (int *)malloc(nlines * sizeof(int));
*amp    = (float *)malloc(nlines * sizeof(float));
*phs   = (float *)malloc(nlines * sizeof(float));
*iq   = (int *)malloc(nlines * sizeof(int));
if(!(*ih)    || !(*ik) || !(*amp)    || !(*phs) || !(*iq) ){
   fprintf(stderr, "\nNo Memory.\n");
   exit(-1);
   }

/* Rewind the file. */
rewind(aphf);

/* Reads the title of the .aph file. */
fgets(line, 200, aphf);

/* Reads frequency components. */
for(indx=0; indx<nlines; indx++) {

    result = fscanf(aphf, "%d %d %f %f %d", 
                          &((*ih)[indx]), &((*ik)[indx]), 
                          &((*amp)[indx]), &((*phs)[indx]), 
                          &((*iq)[indx]));
                          
    fgets(line, 200, aphf); /* Discarding other parts of the line...RMS,CTF? */
    
    if(result==EOF) break;
        
    /* Application of Friedel Inversion. */
    FriedelInversion(&((*ih)[indx]), &((*ik)[indx]), &((*phs)[indx]));
        
    /* Saving the reflection onto the output file. */
/*    printf(" %3d %3d %7.1f %7.1f %2d\n", (*ih)[indx], (*ik)[indx], 
                                         (*amp)[indx], (*phs)[indx], 
                                         (*iq)[indx]);
*/                                         
    }
                                                
fclose(aphf);

return (nlines);
}
/*----------------------------------------------------------------------------*/
void FriedelInversion(int *h, int *k, float *phs)
{
if((*h < 0) || (*h == 0 && *k < 0)){
   *h = - *h;
   *k = - *k;
   *phs = - *phs;
   while(*phs < 0) *phs += 360.0;
   }
}
/*----------------------------------------------------------------------------*/
void Skipline(FILE *f)

/* It skips one line from the input file. */

{
char c = '\0';

while(c!='\n' && !feof(f)) c= fgetc(f);
}
/*----------------------------------------------------------------------------*/
int CompReferencePHS(int Ni, int *hi, int *ki, float *ampi, float *phsi,
                     int *iqi,
                      int Nr, int *hr, int *kr, float *ampr, float *phsr,  
                      int *iqr, float **phsc, float **wgt, float **respot, 
                      float **iptest, float **lspec)
{
int rindx, iindx;
float Af2 = 1.0/(CELL * CELL),
      DstarSQ;
int ncompi,   /* No. comparaciones de una reflexión dada iindx. */
    ncomp;    /* No. total de comparaciones de todas las reflexiones de la
                 imagen de entrada. */

/* Reserva de memoria. */
*phsc   = (float *)malloc(Ni * sizeof(float));
*respot = (float *)malloc(Ni * sizeof(float));
*wgt    = (float *)malloc(Ni * sizeof(float));
*iptest = (float *)malloc(Ni * sizeof(float));
*lspec  = (float *)malloc(Ni * sizeof(float));
if(!(*phsc)    || !(*respot) || !(*wgt)    || !(*iptest) || !(*lspec) ){
   fprintf(stderr, "\nNo Memory.\n");
   exit(-1);
   }


/* Bucle de búsqueda de comparaciones. */
ncomp = 0;
for(iindx=0; iindx<Ni; iindx++){

    /* Resolución. */
    DstarSQ = (hi[iindx] * hi[iindx] * Af2 + ki[iindx] * ki[iindx]* Af2);
    (*respot)[iindx] = 1.0 / sqrt(DstarSQ);
        
    if(DEBUG)
        printf("\n %3d %3d %7.1f %7.1f %2d %7.1f ", 
                   hi[iindx], ki[iindx], ampi[iindx], phsi[iindx], 
                   iqi[iindx], (*respot)[iindx] );

    /* Verificar si se cumplen los límites. */
    if(ampi[iindx] < ZERO || iqi[iindx] > Limits.IQmax ||
       (*respot)[iindx] < Limits.Dresmin || (*respot)[iindx] > Limits.Dresmax) {    
        (*wgt)[iindx] = 0.0;         
        if(DEBUG)
            printf("||                 || "); 
        continue;
        }

    /* Búsqueda de la componente de referencia con la que comparar. */
    ncompi=0;
    for(rindx=0; rindx<Nr && !ncompi; rindx++) {
    
        /* Descartar aquellas cuyo índice (H,K) no es el mismo que el de la
           componente de la imagen de entrada iindx.
           También se descartan aquellas que no cumplen los límites. */
        if(hr[rindx] != hi[iindx]  ||  kr[rindx] != ki[iindx] ||
           ampr[rindx] < ZERO      ||  iqr[rindx] > Limits.IQmax) 
               continue;
        
        /* En este punto se ha encontrado la reflexión equivalente. 
           Podemos actualizar por tanto PHSC. Y a continuación se sale del bucle */
        ncompi++; ncomp++;
        (*phsc)[iindx] = phsr[rindx];
        (*wgt)[iindx] = 1.0;
        (*iptest)[iindx] = 0.0;
        (*lspec)[iindx] = 0.0;
        
        if(DEBUG)
            printf("|| %7.1f %7.1f || %3d %3d %7.1f %2d %7.1f %7.1f", 
                   (*wgt)[iindx], (*phsc)[iindx], hr[rindx], kr[rindx], 
                   ampr[rindx], iqr[rindx], 
                   (*iptest)[iindx], (*lspec)[iindx]);
        
        }
    if(!ncompi && DEBUG)
            printf("||                 || "); 
            
    }                      
return (ncomp);
}    
/*----------------------------------------------------------------------------*/
