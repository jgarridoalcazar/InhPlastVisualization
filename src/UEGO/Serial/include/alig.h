#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/*
 * Now define the data structures, constants, etc for the GHT
 */

#ifndef _ALIGN_
#define _ALIGN_

#define MAXRFL 2500
#define CONSMATH 8
#define DEBUG 0
#define ZERO 0.001
#define CELL 165                /* Dimension of the unit cell: 165 Angstroms. */

/*
 * Data structures
 */
typedef struct {
    float Dresmin, 
          Dresmax;
    int   IQmax;
    } AlignLimits_t;




/*
 * Function prototypes
 */

/*
 * alig Files
 */


/*----------------------------------------------------------------------------*/
float refinephaseorigin(int *,  float, float, 
                        int *, float *, float *,
                        int *, int *, int *, int *,
                        float *, float *, float *, int *, int *,
                        float *);
float phshift(int, int, float, float);
float cmod(float);
/*----------------------------------------------------------------------------*/

/*
 * lectura file
 */


/*----------------------------------------------------------------------------*/
void leerficheros(int filargc, char *filein, char *fileref, int *Ni, int *Ncomp);
void Usage(int, char **);
void Skipline(FILE *);
void FriedelInversion(int *h, int *k, float *phs);
int ReadFile(char *infile, int **ih, int **ik, 
                           float **amp, float **phs, int **iq);
int CompReferencePHS(int Ni, int *hi, int *ki, float *ampi, float *, int *iqi,
                      int Nr, int *hr, int *kr, float *ampr, float *phsr,  
                      int *iqr, float **phsc, float **wgt, float **respot, 
                      int **iptest, int **lspec, AlignLimits_t *);
/*----------------------------------------------------------------------------*/
#endif
