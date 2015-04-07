#ifndef FITNESS_H
#define FITNESS_H

#include <vector>
#include <string>

class Ini;

int MPI_evaluation_loop(Ini * Ini);
double fitness(std::string network_config_file, unsigned int seed, std::vector<std::string> param_names, std::vector<double> param_values);

#endif
