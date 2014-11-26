#include <sys/wait.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#define FITNESS_FILE_NAME "mutual_information"
#define BASE_CONFIG_FILE "../../SimulationConfigTest.cfg"
#define PYTHON_SCRIPT "../LaunchSimulation.py"
#define EXE_FILE "python" // Executed to calculate the fitness



// Simula la red usando el NEST llamado desde python. Configura la red de acuerdo a los valores
// pasados como parámetro a esta función.
// Parametros:
//  simulation_id: string con el nombre de la simulación (se creará un directorio temporal donde se guardará el fichero
// FITNESS_FILE_NAME).
//  param_names: vector<string> con los nombres de los parámetros a utilizar en esta llamada concreta.
//  param_values: vector<double> con los valores de los parámetros a utilizar en la llamada. El número de elementos
// en param_values debe ser el mismo que en param_names.
// Devuelve -2000 si se ha producido un error u otro valor si tuvo éxito
double fitness(std::string simulation_id, std::vector<std::string> param_names, std::vector<double> param_values){

	int status;
	double fitness=-2000;

	std::cout << "Running fitness function" << std::endl;
	std::string command_param = std::string(EXE_FILE) + std::string(" ") + std::string(PYTHON_SCRIPT);
	command_param += std::string(" -c ") + std::string(BASE_CONFIG_FILE);
	command_param += std::string(" simulation.simulation_name=") + simulation_id;
	for (unsigned int i=0; i<param_names.size(); ++i){
		std::ostringstream s;
		s << param_names[i] << "=" << param_values[i];
		command_param += " " + s.str();
	}

	std::cout << "Running the simulations: " << command_param << std::endl;
	// Ejecutar la simulación
	system(command_param.c_str()); // Launch a subprocess to run the simulation

    // Leer el fichero con el valor de la función fitness
    std::ifstream file_id;
    std::string file_name = "./results/" + simulation_id + "/" + FITNESS_FILE_NAME;
    file_id.open(file_name.c_str());
    file_id >> fitness;
    file_id.close();
    if (!file_id){
    	perror("Reading fitness");
    }

	printf(">returned fitness: %lf\n",fitness);
	return(fitness);

}

