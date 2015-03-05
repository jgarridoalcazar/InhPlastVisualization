#ifndef UEGO_H
#define UEGO_H

////////////////////////////////////////////////////////////
// $Id: uego.h,v 2.6 1998/04/01 21:05:18 jelasity Exp $
// uego.h
// collects different headers for simplification of usage;
// should be included from all modules
////////////////////////////////////////////////////////////
// modification history:
//	jelasity 98 02 20 ini.h -> uegoini.h
//		98 03 10 nreal.h nbin.h
////////////////////////////////////////////////////////////

#define UEGO_VERSION "2.6"

// --- for sending errors and warnings, must be defined in a user interface
void message( char*, short );

// --- defines for message
#define MSG_INFORMATION  2
#define MSG_ERROR        1
#define MSG_NOTHING      0


// --- define to keep source clean
#define INI (Master::ini())
#include <mpi.h>
#include <stdio.h>
#include <math.h>

class SearchSpElement;
class SpeciesList;
class Master;

#include "uegoini.h"
#include "uegorand.h"
#include "master.h"

#include "searchsp.h"
#include "nreal.h"
#include "nbin.h"
#include <pthread.h>
#include "speclist.h"

#endif
