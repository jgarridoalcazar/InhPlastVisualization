////////////////////////////////////////////////////////////
// $Id: fitness.cc,v 2.0 24/03/2015 $
// fitness.cc
// definition of slave nodes functioning
// Created by Jesus Garrido.
////////////////////////////////////////////////////////////
//------------------------------------------------------------------

#include <sys/wait.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <Python.h>
#include <math.h>
#include <uego.h>

#define PYTHON_FILE_NAME "./src/CInterface.py"
#define MODULE_NAME "CInterface"
#define FUNCTION_NAME "RunSimulation"

// Receive parameter configurations through MPI interface,
// evaluate the fitness function and send the evaluation result
// back.
int MPI_evaluation_loop(Ini * Ini){

	bool finished = false;

	// Reception buffer NDIM*double + unsigned int (seed)
	int Nbytes = Ini->Dimension()*sizeof(double) + sizeof(unsigned int);
	void * recvBuf = (void *) new char [Nbytes];
	double * x_norm_values = new double [Ini->Dimension()];
	unsigned int seed;
	int position;
	double sim_result;

	// Get the param names and values vectors
	std::vector<std::string> param_names(Ini->Dimension());
	std::vector<double> param_values(Ini->Dimension());

	MPI_Status	 status;

	while (!finished){
		MPI_Recv(recvBuf, Nbytes, MPI_PACKED, 0, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
		// Check the tag of the message to now whether the simulation is finished
		if (status.MPI_TAG==1){
			char msg [150];
			sprintf( msg, "WORKER :: Received ending signal from the master");
			message( msg, MSG_INFORMATION );
			finished = true;
			continue;
		}


		position = 0;
		// Unpack the seed and the normalized values
		MPI_Unpack(recvBuf,Nbytes,&position,x_norm_values,Ini->Dimension(),MPI_DOUBLE,MPI_COMM_WORLD);
		MPI_Unpack(recvBuf,Nbytes,&position,&seed,1,MPI_UNSIGNED,MPI_COMM_WORLD);

		//---Allocation memory
		for(int i=0;i<Ini->Dimension();i++){
			param_names[i]=Ini->ParameterName(i);
			if (Ini->ParameterScale(i)==LOGARITHMIC){
				double logmin = log10(abs(Ini->Lowb(i)));
				double logmax = log10(abs(Ini->Upb(i)));
				param_values[i]= pow(10.0,(x_norm_values[i]*(logmax - logmin)))*Ini->Lowb(i);
			} else {
				param_values[i]=x_norm_values[i]*(Ini->Upb(i)-Ini->Lowb(i)) + Ini->Lowb(i);
			}
		}

		// Evaluate the fitness function
		sim_result = fitness(std::string(Ini->NetworkConfigFile()),seed,param_names, param_values);

		// Send the fitness function result back to the master
		if(MPI_Send(&sim_result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD) != MPI_SUCCESS){
			printf("An error in MPI_Send() Slave :: send simulation result to master\n");
			finished = true;
		}

	}

	delete [] (char *) recvBuf;
	delete [] x_norm_values;


	return 0;
}




// Simula la red usando el NEST llamado desde python. Configura la red de acuerdo a los valores
// pasados como parámetro a esta función.
// Parametros:
//	network_config_file: File name of the network configuration file
//  seed: The seed used to initialize the simulation.
//  param_names: vector<string> con los nombres de los parámetros a utilizar en esta llamada concreta.
//  param_values: vector<double> con los valores de los parámetros a utilizar en la llamada. El número de elementos
// en param_values debe ser el mismo que en param_names.
// Devuelve -2000 si se ha producido un error u otro valor si tuvo éxito
double fitness(std::string network_config_file, unsigned int seed, std::vector<std::string> param_names, std::vector<double> param_values){

	PyObject *py_imp_str, *py_imp_handle, *py_imp_dict, *py_imp_load_source, *py_lib_mod, *py_args_tuple, *pValue, *py_key_name, *py_param_value; // New references
	PyObject *py_dir, *py_lib_name, *py_lib_func, *pArgs; //stolen
	PyObject *py_lib_mod_dict, *py_func; // Borrowed references
	int i;
	double result = 0.0;

	char msg[1000];
	sprintf(msg, "Simulating with seed %ld",seed);
	for (i=0; i<param_names.size(); ++i){
		sprintf(msg, "%s, %s=%e", msg, param_names[i].c_str(), param_values[i]);
	}
	sprintf(msg, "%s\n", msg);
	message(msg, MSG_INFORMATION);

	//import our python script using the imp library (the normal import doesn't allow you to grab py files in random directories)
	py_imp_str = PyString_FromString("imp");
	py_imp_handle = PyImport_Import(py_imp_str); //normal python import for imp
	Py_DECREF(py_imp_str);

	if (py_imp_handle != NULL) {
		py_imp_load_source = PyObject_GetAttrString(py_imp_handle, "load_source");
		// py_imp_load_source is a new reference */

		if (py_imp_load_source && PyCallable_Check(py_imp_load_source)) {

			py_dir = PyString_FromString(PYTHON_FILE_NAME);
			if (!py_dir) {
				Py_DECREF(py_imp_handle);
	            fprintf(stderr, "Cannot convert string argument (directory - %s)\n",PYTHON_FILE_NAME);
	            return 0;
	        }

	    	py_lib_name = PyString_FromString(MODULE_NAME);
			if (!py_lib_name) {
				Py_DECREF(py_imp_handle);
	            fprintf(stderr, "Cannot convert string argument (module - %s)\n",MODULE_NAME);
	            return 0;
	        }

	        py_lib_func = PyString_FromString(FUNCTION_NAME);
			if (!py_lib_func) {
				Py_DECREF(py_imp_handle);
	            fprintf(stderr, "Cannot convert string argument (function - %s)\n",FUNCTION_NAME);
	            return 0;
	        }

	    	//setup args for imp.load_source
			py_args_tuple = PyTuple_New(2);
	    	PyTuple_SetItem(py_args_tuple, 0, py_lib_name); //stolen
			PyTuple_SetItem(py_args_tuple, 1, py_dir); //stolen

		    //call imp.load_source
			py_lib_mod = PyObject_CallObject(py_imp_load_source, py_args_tuple);
			Py_DECREF(py_args_tuple);
	        if (py_lib_mod == NULL) {
	        	Py_DECREF(py_imp_load_source);
	            Py_DECREF(py_imp_handle);
	            PyErr_Print();
	            fprintf(stderr,"Call to load_source failed\n");
	            return 0;
	        }

			//get function object
			py_lib_mod_dict = PyModule_GetDict(py_lib_mod); //borrowed

			//get function object
			py_func = PyDict_GetItem(py_lib_mod_dict, py_lib_func);

			if (py_func && PyCallable_Check(py_func)) {
				py_args_tuple = PyTuple_New(1);
			    pArgs = PyDict_New();

			    // Insert configuration file name
			    py_key_name = PyString_FromString("config_file_name");
				if (!py_key_name) {
					Py_DECREF(pArgs);
	    			Py_DECREF(py_imp_handle);
	    			Py_DECREF(py_lib_mod);
		            fprintf(stderr, "Cannot convert config_file_name string\n");
	    	        return 0;
	            }

	            py_param_value = PyString_FromString(network_config_file.c_str());
				if (!py_param_value) {
					Py_DECREF(pArgs);
	    			Py_DECREF(py_imp_handle);
	    			Py_DECREF(py_lib_mod);
		            fprintf(stderr, "Cannot convert config file name string %s\n",network_config_file.c_str());
		            return 0;
	            }

	            PyDict_SetItem(pArgs, py_key_name, py_param_value);
			    Py_DECREF(py_key_name);
			    Py_DECREF(py_param_value);

			    // Insert int simulation seed
			    py_key_name = PyString_FromString("simulation.seed");
				if (!py_key_name) {
					Py_DECREF(pArgs);
	    			Py_DECREF(py_imp_handle);
	    			Py_DECREF(py_lib_mod);
		            fprintf(stderr, "Cannot convert simulation.seed string\n");
	    	        return 0;
	            }

	            py_param_value = PyInt_FromLong(seed);
				if (!py_param_value) {
					Py_DECREF(pArgs);
	    			Py_DECREF(py_imp_handle);
	    			Py_DECREF(py_lib_mod);
		            fprintf(stderr, "Cannot convert int parameter\n");
		            return 0;
	            }

	            PyDict_SetItem(pArgs, py_key_name, py_param_value);
			    Py_DECREF(py_key_name);
			    Py_DECREF(py_param_value);

			    for (i=0; i<param_names.size();i++){
			    		// Insert float param value
			    		py_key_name = PyString_FromString(param_names[i].c_str());
			    		if (!py_key_name) {
			    			Py_DECREF(pArgs);
			    			Py_DECREF(py_imp_handle);
			    			Py_DECREF(py_lib_mod);
			    			fprintf(stderr, "Cannot convert mfgocsynapsis.minus_plus_ratio string\n");
			    			return 0;
			    		}

			    		py_param_value = PyFloat_FromDouble(param_values[i]);
			    		if (!py_param_value) {
			    			Py_DECREF(pArgs);
			    			Py_DECREF(py_imp_handle);
			    			Py_DECREF(py_lib_mod);
			    			fprintf(stderr, "Cannot convert float parameter\n");
			    			return 0;
			    		}

			    		PyDict_SetItem(pArgs, py_key_name, py_param_value);
						Py_DECREF(py_key_name);
						Py_DECREF(py_param_value);

			    }



			    PyTuple_SetItem(py_args_tuple, 0, pArgs); // pArgs stolen

			    pValue = PyObject_CallObject(py_func, py_args_tuple);
			    //pValue = (PyObject *) 1;
	            Py_DECREF(py_args_tuple);
	            if (pValue != NULL) {
	            	//result = 0.0;
	            	result = PyFloat_AsDouble(pValue);
	            	Py_DECREF(pValue);
	            }else {
	                Py_DECREF(py_imp_handle);
	    			Py_DECREF(py_lib_mod);
	                PyErr_Print();
	                fprintf(stderr,"Call failed\n");
	                return 0;
	            }
	        }else{
	        	if (PyErr_Occurred()){
	        		PyErr_Print();
	        	}
	        	fprintf(stderr, "Cannot find function \"%s\"\n", FUNCTION_NAME);
	        }

	        Py_DECREF(py_lib_mod);
	    }else{
	    	if (PyErr_Occurred()){
	    		PyErr_Print();
	    	}
	        fprintf(stderr, "Cannot find function \"%s\"\n", "load_source");
	    }

	    Py_DECREF(py_imp_load_source);
	    Py_DECREF(py_imp_handle);
	}else {
		PyErr_Print();
	    fprintf(stderr, "Failed to load \"%s\"\n", "imp");
	    return 0;
	}

	sprintf(msg, "Simulation ended with seed %ld",seed);
	for (i=0; i<param_names.size(); ++i){
		sprintf(msg, "%s, %s=%e", msg, param_names[i].c_str(), param_values[i]);
	}
	sprintf(msg, "%s. Returned: %f\n", msg, result);
	message(msg, MSG_INFORMATION);

	return result;
}

