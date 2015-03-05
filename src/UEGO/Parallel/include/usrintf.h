#ifndef USRINTF_H
#define USRINTF_H

////////////////////////////////////////////////////////////
// $Id: usrintf.h,v 2.7 1998/04/01 21:05:18 jelasity Exp $
// usrintf.h
// defines for the main modul of the user interface
////////////////////////////////////////////////////////////
// modification history:
//	jelasity 98 01 17 DEFAULT_MESSAGE, VERSION
//	jelasity 98 01 24 DEFAULT_MESSAGE, VERSION
////////////////////////////////////////////////////////////


#define UEGO_ININAME "uegoini"
#define DEFAULT_MESSAGE "\
uego <ext> [-T<filebase>] [-r<repeatcount>] [-s]\n\
	Optimalization using the ini file 'uegoini.<ext>'. Saves to stdout.\n\
	-T tracing; e.g. '-T./d' means that the checkpoint files will be\n\
		placed to ./d0000.ckp, ./d0001.ckp and so on.\n\
	-r repeat the experiment. e.g. '-r5' results in 5 independent\n\
		experiments.\n\
	-s save result in short format (no species, only summary)\n\n\
uego -c<ext>\n\
	setup uego interactively (creating or updating uegini.<ext>)\n\n\
uego -c<ext> [-s<>] [-p<>] [-N<>] [-M<>] [-t<>] [-l<>] [-r<>] [-x<>]\n\
	setup uego from command line (creating or updating uegini.<ext>)\n\
	Parameter order is important. The parameters:\n\
	-s: random seed, -p: parameter file, -N: max. number of evaluations,\n\
	-M: max. number of species, -t: treshold, -l: levels,\n\
	-r: smallest radius, -x: parameter to ignore\n"


#endif

