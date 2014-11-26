#ifndef UEGOCONF_H
#define UEGOCONF_H

////////////////////////////////////////////////////////////
// $Id: uegoconf.h,v 2.6 1998/04/01 21:05:18 jelasity Exp $
// uegoconf.h
// contains the defines for uego configuration (Ini::Ini())
// included from ini.cc and configur.cc only
////////////////////////////////////////////////////////////
// modification history:
//
////////////////////////////////////////////////////////////


// ----- handling search space types --------------------------------------


#define NUM_OF_TYPES 2

#define	NDIMREAL_TYPE 0
#define	BINARY_TYPE 1

#define VALID_TYPE(x) ( (x)==NDIMREAL_TYPE || (x)==BINARY_TYPE )
#define BOUND_NEEDED(x) ( (x)==NDIMREAL_TYPE )

#define TYPE_HE "\n\
Currently, 2 types are supported:\n\
\tn dimensional real space: 0\n\
\tn dimensional binary space: 1\n\
give one of the above numbers.\n\n"


// ----- tokens for values ---------------------------------------------
// one token should not include another!

// --- random numbers
#define SEED_T "random seed:"

// --- objective function
#define TYPE_T "search space type:"
#define FNUM_T "function index:"
#define REC_P "number of entering firms:"
#define DIM_T "dimension of space:"
#define PARNUM_T "number of constant parameters:"
#define PAR_T "constant parameter values:"
#define LOWB_T "lower bound vector:"
#define UPB_T "upper bound vector:"

// --- uego special
#define MAXEV_T "function evaluations:"
#define MAXSPEC_T "max. species number:"
#define THR_T "stability threshold:"
#define LEVELS_T "levels:"
#define LAST_R_T "smallest radius:" // for Configure only
// --- automatic
#define NEWSPEV_T "function evaluations in new species:"
#define EVALS_T "evaluations for levels:"
#define RADII_T "radii for levels:"

// --- Parameters file
#define PARFILE_T "parameters file: "
#define NUM_DEMAND_T "num_of_demand_points(n)"
#define NUM_EXIT_FAC_T "num_of_existing_facilities(m)"
#define LENGHT_CHAIN_T "length_of_chain(k)"
#define BOUNDS_VAR_T "Bounds_for_the_variables_x1_x2_alpha:"

#define PARAMETER_G_FUNCTION "Parameter_for_g_function(lambda__or__beta__in_hodson):"
#define A0 "a0"
#define A1 "a1"
#define B1 "b_1"
#define B2 "b_2"
#define P "p_"
#define PARAM_F "Parameter_for_F:"
#define PARAMETERS_PSI "Parameters_for_Psi_functions:"
#define LOC_EXIT_FAC "location_of_the_existing_facilities(f[m][2]):"
#define LOC_DEMAND_POINT "location_of_the_demand_points(p[n][2]):"

#define BUYING_POWER "buying_power_at_the_demand_points(w[n]):"
#define WEIGHT_OF_ATTRAC "weight_of_attractiveness_to_x_of_the_demand_points(v[n]--gamma):"
#define ATTRAC_EXIST_FAC "attractiveness_of_demand_points_to_the_existing_facilities(a[n][m]--alfa):"

// ----- help messages ---------------------------------------------

#define STARTUP_HE "\n\
Welcome to uego configuration!\n\n\
Please type in the requested values. If the default value in () is ok\n\
just type <enter>. If you are editing an existing configuration, then\n\
the defalult values are read from that configuration.\n\
If not sure, type h<enter>; a help message will be shown reagarding the\n\
actual item if available. If you want to cancel configuration, type q<enter>.\n\n"

#define SEED_HE "\n\
The random seed determines the random numbers created during optimalization.\n\
With the same random seed, the same results will come out on every platform.\n\
Choose a number from [ 1 , 2^32-1 ], large numbers are better.\n\n"

//#define TYPE_HE : see above at search spaces

#define FNUM_HE "\n\
It is possible to encode several functions so that you could run experiments\n\
on them without recompiling the system over and over again. The number you\n\
give here can be seen from the evaluator function, and can be used in e.g.\n\
a 'switch' command. If you did not use this feature, this value is not\n\
effective; otherwise give the number of function you would like to experiment\n\
with.\n\n"

#define REC_HE "\n\
Give the number of new facilities you want to localize\n\n"


#define DIM_HE "\n\
Give the dimensionality of the search space you are using. In the real case,\n\
it is the number of real parameters of the function, in the binary case, the\n\
number of bits, etc.\n\
It must be at least 1.\n\n"

#define PARS_HE "\n\
It is possible to pass constant parameters to a function. Here, you can\n\
give the name of a text file to load these parameters from. Hitting <enter>\n\
will cancel choosing a file.\n\
The file format is simple: it should contain the number of parameters and\n\
then the real parameters separated by whitespace. Lines beginning with a '#'\n\
are ignored. Lines should be no longer than 10000 characters. Be sure to\n\
close the last line with a <newline> character.\n\n"

#define PARNUM_HE "\n\
It is possible to pass constant parameters to a function. Here, you can give\n\
the number of these parameters. If your function does not use any, give 0.\n\n"

// --- uego special
#define MAXEV_HE "\n\
Give the number of evaluations you want to devote to this experiment.\n\
The actual evaluation number will not exceed the given limit but may\n\
be less.\n\n"

#define MAXSPEC_HE "\n\
The maximum number of species allowed to exist in the same time.\n\
This value can be understood as a kind of exploration parameter: greater\n\
values mean more extensive exploration. If your function does not have\n\
too many 'peaks' give small values, otherwise give a larger value.\n\
It is safe to give smaller values than the actual number of local optima\n\
but if a large value is given for a function with few local optima then most\n\
of the evaluations will be devoted to finding the non-existing peaks.\n\
Do not forget, that the number of function evaluations is shared among the\n\
existing species, so this value should be much smaller than the number of\n\
function evaluations, of course.\n\n"

#define THR_HE "\n\
This is a threshold that represents a trade-off between resolution (or\n\
exploration) and stability. Setting 1 means that every species at every\n\
level is guaranteed to receive enough function evaluations to move from its\n\
creation place as far as the diameter of the space. 0.5 means the half of\n\
this distance and so on.\n\
Setting 1 seems to be a good choice in general.\n\n"

#define LEVELS_HE "\n\
Quite hard to explain... Well if your function has a kind of fractal-like\n\
structure, it may be reasonable to set a higher value; the role of this\n\
parameter is still under investegation (experiments show that the performance\n\
may depend on this value, sometimes higher values are slighly better, but\n\
sometimes not). Setting e.g. 3 seems to be safe in general.\n\n"

#define LAST_R_HE "\n\
The resolution of search. The radius of species on the final level.\n\
The optima found are at least this far from each other. The smaller the\n\
given value is the more evaluaions are needed.\n\
This parameter is NOT directly related to the expected error of the\n\
solutions; set very small values only if you want to find different local\n\
optima that are close to each other, or if you are looking for needle-like\n\
peaks.\n\n"

#define MISS_HE "\n\
Every 1 of the last 5 parameters is determined by the remaining 4. Here\n\
you can specify the one that you want to be counted using the other 4\n\
values. It is reasonable to choose the most obscure one if you do not have\n\
a better idea. Currently, the valid choises are:\n\n\
   0: function evaluations\n\
   1: max. species number\n\
   2: stability threshold\n\
   4: smallest radius\n\n\
If you got an error message after this (something like 'No solution exists')\n\
then try one (or more) of these next time: increase function evaluations,\n\
decrease threshold, increase smallest radius, decrease max species number.\n\n"


// --- ini out header ------------------------------------------------------

#define INI_OUT_HEADER "\
#\n\
# This is an automatically generated UEGO configuration file.\n\
# Edit only if you really know, what you're doing!\n\
#\n# For configuration, use 'uego -c' instead.\n\
#\n"

#endif
