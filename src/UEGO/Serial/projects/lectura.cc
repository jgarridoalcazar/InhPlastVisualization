#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "uego.h"
#include "alig.h"

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
      printf("\n error abriendo el fichero\n");                          /* Unable to open input file. */
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
/*----------------------------------------------------------------------------*/
int CompReferencePHS(int Ni, int *hi, int *ki, float *ampi, float *phsi,
                     int *iqi,
                      int Nr, int *hr, int *kr, float *ampr, float *phsr,  
                      int *iqr, float **phsc, float **wgt, float **respot, 
                      int **iptest, int **lspec, AlignLimits_t *Limits)
{
int rindx, iindx;
float Af2 = 1.0/(CELL * CELL),
      DstarSQ;
int ncompi,   /* No. comparaciones de una reflexin dada iindx. */
    ncomp;    /* No. total de comparaciones de todas las reflexiones de la
                 imagen de entrada. */

/* Reserva de memoria. */
*phsc   = (float *)malloc(Ni * sizeof(float));
*respot = (float *)malloc(Ni * sizeof(float));
*wgt    = (float *)malloc(Ni * sizeof(float));
*iptest = (int *)malloc(Ni * sizeof(int));
*lspec  = (int *)malloc(Ni * sizeof(int));
if(!(*phsc)    || !(*respot) || !(*wgt)    || !(*iptest) || !(*lspec) ){
   fprintf(stderr, "\nNo Memory.\n");
   exit(-1);
   }


/* Bucle de bsqueda de comparaciones. */
ncomp = 0;
for(iindx=0; iindx<Ni; iindx++){

    /* Resolucin. */
    DstarSQ = (hi[iindx] * hi[iindx] * Af2 + ki[iindx] * ki[iindx]* Af2);
    (*respot)[iindx] = 1.0 / sqrt(DstarSQ);
    (*wgt)[iindx] = 0.0;
    (*iptest)[iindx] = 0;
    (*lspec)[iindx] = 0;
        
    if(DEBUG)
        printf("\n %3d %3d %7.1f %7.1f %2d %7.1f ", 
                   hi[iindx], ki[iindx], ampi[iindx], phsi[iindx], 
                   iqi[iindx], (*respot)[iindx] );

    /* Verificar si se cumplen los l?ites. */
    if(ampi[iindx] < ZERO || iqi[iindx] > Limits->IQmax ||
       (*respot)[iindx] < Limits->Dresmin || (*respot)[iindx] > Limits->Dresmax) {    
        (*wgt)[iindx] = 0.0;         
        if(DEBUG)
            printf("||                 || "); 
        continue;
        }

    /* Bsqueda de la componente de referencia con la que comparar. */
    ncompi=0;
    for(rindx=0; rindx<Nr && !ncompi; rindx++) {
    
        /* Descartar aquellas cuyo ?dice (H,K) no es el mismo que el de la
           componente de la imagen de entrada iindx.
           Tambi? se descartan aquellas que no cumplen los l?ites. */
        if(hr[rindx] != hi[iindx]  ||  kr[rindx] != ki[iindx] ||
           ampr[rindx] < ZERO      ||  iqr[rindx] > Limits->IQmax) 
               continue;
        
        /* En este punto se ha encontrado la reflexin equivalente. 
           Podemos actualizar por tanto PHSC. Y a continuacin se sale del bucle */
        ncompi++; ncomp++;
        (*phsc)[iindx] = phsr[rindx];
        (*wgt)[iindx] = 1.0;
        (*iptest)[iindx] = 0;
        (*lspec)[iindx] = 0;
        
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
