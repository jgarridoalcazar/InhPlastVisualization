#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXRFL 2500
#define CONSMATH 8

/*----------------------------------------------------------------------------*/
void refinephaseorigin_(int *, float *, int *, int *, float *, float *, float *,
                        int *, int *, int *, int *, int *, float *, float *,
                        int *, float *, float *, int *, int *, int *,
                        float *, float *, float *, int *, int *,
                        float *, float *, float *, short int [][121]);
float phshift(int, int, float, float, float, float, float *, float);
float cmod(float);
/*----------------------------------------------------------------------------*/
void refinephaseorigin_(
    int   *iboxphs,
    float *step,
    int   *iorigt,
    int   *in1,
    float x0,
    float x1,
    float *errmin,
    int   *ishmin,
    int   *iskmin, 
    int   *ninec, 
    int   *ncomp,
    int   *ncompi,
    float *dresmin,
    float *dresmax,
    int   *iqmax,
    float *tilth,
    float *tiltk,
    int   *iih,
    int   *iik,
    int   *iqin,
    float *phsi,
    float *phsc,
    float *wgt,
    int   *iptest,
    int   *lspec,
    float *respot,
    float *beamshft,
    float *bsh,
    short int irp[121][121]
    )
{
int ncompt, ncomp9,
    ish, isk, irefc;
float sh, sk, sk0, serr, p, delta, serr9;
int iserr;


printf("   ORIGIN REFINEMENT DONE BETWEEN %5d OF THE NEW REFLECTIONS\n"
       "  AND ,%5d OF THE REFLECTIONS FROM PREVIOUS FILMS. \n"
       "   %5d REFLECTIONS HAVE PHASES CONSTRAINED BY SYMMETRY.\n",
       *ncompi, *ncomp, *ninec);
ncompt = *ninec + *ncomp;
ncomp9 = *ninec + *ncompi;
if(*iorigt == 1) 
    printf(" FOR \"EQUAL WEIGHT PER INPUT SPOT ON NEW FILM\" TYPE\n"
           "   OF ORIGIN REFINEMENT, TOTAL COMPARISONS =%5d\n", ncomp9);
printf(" THIS IS A TOTAL OF %5d COMPARISONS.\n", ncompt);


*errmin=MAXFLOAT;
 serr  = 0.0;    

/*-------------------------------------------       
// Relacion entre sh, sk y x[0], x[1]
//sh = sk0 = -0.5 * (*iboxphs + 1) * (*step);
 // for (ish=0; ish<*iboxphs; ish++){
//    sh += *step;
//    sk = sk0;
//     for(isk=0; isk<*iboxphs; isk++){
//        sk += *step;
----------------------------------------------*/
                            /* Origtilt Error Funcion. */
sh=x0;
sk=x1;
        for(irefc=0; irefc<*in1; irefc++){
        
            /* Discarding reflection. */
            if(iqin[irefc] > *iqmax ||
                respot[irefc] < *dresmin ||
                respot[irefc] > *dresmax) continue;
            
             p = phsi[irefc] + phshift(iih[irefc],iik[irefc],sh,sk);
                                      
            /*------------------------  Origtilt ---------------------------- */
            if(wgt[irefc] >= 1.0){
                delta = cmod(p - phsc[irefc]);
                if(delta < 0.0) delta = -delta;
                serr += wgt[irefc] * delta;        /* Origtilt Error Funcion. */
                }
            if(!lspec[irefc] || 
                iqin[irefc] > *iqmax ||
                respot[irefc] < *dresmin ||
                respot[irefc] > *dresmax) continue;
            p = cmod(p);
            if(p < 0.0) p = -p;
            delta = p - iptest[irefc];
            if(delta < 0.0) delta = -delta;
            if(delta > 90.0) delta = 180.0 - delta;
            serr += delta;
            /*----------------------------------------------------------------*/
            }
        if(ncomp9) serr9 = serr / ncomp9;
        if(ncompt) serr /= ncompt;
        if(*iorigt == 1) serr  = serr9;
        if(serr < *errmin){
 /*-----------------------------------------
            *skmin = sk;
            *shmin = sh;

            *iskmin = isk + 1;
            *ishmin = ish + 1;
 ---------------------------------------------*/
            *errmin = serr;
            }
        iserr = (int)(90.0 - serr)/10;
        if(iserr < 0) iserr = 0;
  //      irp[isk][ish] = iserr;               /* indices have been swapped. */
        }
    }
    
return;
}
/*----------------------------------------------------------------------------*/
float phshift(int ih, int ik, float ox, float oy)
{
#define astar beamshft[1]
#define bstar beamshft[2]

return(ih*ox + ik*oy);
}
/*----------------------------------------------------------------------------*/
float cmod(float p)
/* Returns p module 360 in [-180,180] */
{
p = fmodf(p,360.0); 
if(p > 180.0) p -= 360.0;
if(p <= -180.0) p += 360.0;
if(p <= -179.9) p = 180.0;
return(p);
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

