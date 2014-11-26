#include "uego.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fitness.c"
#include <vector>

// INI.Lowb(i) and INI.Upb(i) are lower and upper bounds on x[i] max(i)=dim-1;
double	NDimRealElement::Value()
{
std::vector<std::string> param_names(3);
param_names[0] = "gocgocsynapsis.learning_step";
param_names[1] = "gocgocsynapsis.minus_plus_ratio";
param_names[2] = "gocgocsynapsis.max_weight";

std::vector<double> param_values(3);
param_values[0] = 5.0e-3;
param_values[1] = 0.5;
param_values[2] = -1.0e-8;

return(fitness("simulation_test",param_names, param_values));
getchar();

}

