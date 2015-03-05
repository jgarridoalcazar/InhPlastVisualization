#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "uego.h"
#include "alig.h"



float refinephaseorigin(
    int   *in1,
    float x0,
    float x1,
    int   *ncomp,
    float *dresmin,
    float *dresmax,
    int   *iqmax,
    int   *iih,
    int   *iik,
    int   *iqin,
    float *phsi,
    float *phsc,
    float *wgt,
    int   *iptest,
    int   *lspec,
    float *respot
    )
{
int   irefc;
float sh, sk, serr, p, delta;

 serr  = 0.0;    

                            /* Origtilt Error Funcion. */
sh=x0;
sk=x1;
 
        for(irefc=0; irefc<*in1; irefc++){

            /* Discarding reflection. */
            if(iqin[irefc] > *iqmax ||
                respot[irefc] < *dresmin ||
                respot[irefc] > *dresmax) continue;

             p = phsi[irefc] + phshift(iih[irefc],iik[irefc],sh,sk);

            /*------------------------  Origtilt --------------------------- */
            if(wgt[irefc] >= 1.0){
                delta = cmod(p - phsc[irefc]);
                if(delta < 0.0) delta = -delta;
                serr += wgt[irefc] * delta;       /* Origtilt Error Funcion. */
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
     serr /= *ncomp;
    
return serr;
}
/*----------------------------------------------------------------------------*/
float phshift(int ih, int ik, float ox, float oy)
{
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

