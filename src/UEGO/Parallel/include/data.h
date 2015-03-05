#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
 * Now define the data structures, constants, etc for the GHT
 */

#ifndef _GHT_H_
#define _GHT_H_

/*
 * Constants
 */

/* ImageData constants */

#define Bdata 0
#define Idata 1
#define Fdata 2

/* TableData constants */

#define Model 0
#define Image 1

/* PGM images constants */

#define PGM_MAGIC1 'P'
#define PGM_MAGIC2 '2'
#define RPGM_MAGIC2 '5'
#define PGM_FORMAT (PGM_MAGIC1 * 256 + PGM_MAGIC2)
#define RPGM_FORMAT (PGM_MAGIC1 * 256 + RPGM_MAGIC2)

/* GHT constants */

#define MAX_LENGTH 20 /* Factor de reasignacion, revisar realloc */
#define MAX_DISTANCES 350 /* Distancia maxima entre los puntos de una pareja */
#define MAX_SCALES 200 /* Numero de escalas distintas que se exploraran */

/* Mathematics constants */

#if defined(__STDC__) || defined(ANSI)

#define PI 3.1415927f
#define Rad2Deg 57.29577866f

#else /* Traditional K&R */

#define PI 3.1415927
#define Rad2Deg 57.29577866

#endif

/*
 * Data structures
 */

typedef unsigned char Byte;

typedef struct {
	short width, height;
	int type;
	Byte *bdata;
	int *idata;
	float *fdata;
} ImageData;

typedef struct {
	short width, height;
	short num_points;
	short *coorx;
	short *coory;
	float *gradx;
	float *grady;
	float *angle;
} ContourPoints;

typedef struct {
	short index;
	short length;
	int *point1;
	int *point2;
} PairPoints;

typedef struct {
	short width, height;
	int type;
	int *idata;
	int **distances;
	PairPoints *pairs;
} TableData;

/*
 * Function prototypes
 */

/*
 * PGM files
 */

ImageData *ReadPGMimage(FILE *fp);
void WritePGMimage(FILE *fp, ImageData *image);

/*
 * Data management
 */

ImageData *AllocID(short type, /* data type: byte -> Bdata, int -> Idata, float -> Fdata */
				   short width, short height,	/* dimensions */
				   Byte b_init, int i_init, float f_init /* init values for the data */
				   );
ContourPoints *AllocCP(short num_points, /* Number of points in the contour */
					   short width, short height /* dimensions of the image */
					   );
TableData *AllocTD(short width, short height /*dimensions of the table */
				   );
void FreeID(ImageData *ID);
void FreeCP(ContourPoints *CP);
void FreeTD(TableData *TD);
ImageData *CopyID(ImageData *ID);
ContourPoints *CopyCP(ContourPoints *CP);
TableData *CopyTD(TableData *TD);
ContourPoints *CPconcat(ContourPoints *cp1, ContourPoints *cp2);
ImageData *MergeID(ImageData *id1, ImageData *id2);
ImageData *ID2IDbyte(ImageData *ID);
ImageData *CP2IDbyte(ContourPoints *CP, Byte value);
ImageData *TD2IDbyte(TableData *TD);
ContourPoints *ID2CP(ImageData *edges, ImageData *gradx, ImageData *grady);

/*
 * Transformation routines
 */

ContourPoints *HomoCP(ContourPoints *input, int despx, int despy, float rotation, float scale);
ContourPoints *DeformCP(ContourPoints *input, float **xi_x, float **xi_y, ContourPoints *deform);

/*
 * Edge computation
 */

ContourPoints *Canny(ImageData *input, float sigma, int threshold);

/*
 * GHT routines
 */

TableData *ComputeTD(ContourPoints *edges, int resolution, int chi);
int WNCCrotation(TableData *templ, TableData *image, int window, float *accuracy);
int WNCCscale(TableData *templ, TableData *image, int rotation, int window);
float ExtractPoints(ContourPoints *edges1, TableData *table1, ContourPoints *edges2, TableData *table2, int rotation, int scale);
ImageData *ExtractCenter(ContourPoints *edges1, TableData *table1, ContourPoints *edges2, TableData *table2, int rotation, int scale, float *accuracy);

/*
 * Useful routines
 */

float ComputeAngle(float x, float y);
float gaussian(float x, float s);
float hypotenuse(float x, float y);

#endif
