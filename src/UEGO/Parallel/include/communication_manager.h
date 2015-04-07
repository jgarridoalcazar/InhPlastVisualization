#ifndef COMMUNICATION_MANAGER_H
#define COMMUNICATION_MANAGER_H

#include "thread.h"
#include "wqueue.h"
#include <list>
#include <vector>
#include <utility>	// std::pair
#include <map>
#include <time.h>
#include "uego.h"

class SearchSpElement;
class Master;

class CommunicationManager : public Thread
{
    wqueue<SearchSpElement *>& input_queue;
    wqueue<SearchSpElement *>& output_queue;
    bool & end_simulation;
 
  	public:
    	CommunicationManager(wqueue<SearchSpElement *>& to_simulate_queue, wqueue<SearchSpElement *>& finished_queue, bool& simulation_end) :
    		input_queue(to_simulate_queue), output_queue(finished_queue), end_simulation(simulation_end) {};
 
    void* run() {

    	MPI_Request request;
    	MPI_Status status;

    	end_simulation = false;

    	int simulation_finished = 0;

    	unsigned int endedProcesses = 1;

    	double result;

    	unsigned int Nbytes = INI.Dimension()*sizeof(double) + sizeof(unsigned int);
    	void * sendBuf = new char [Nbytes];
    	double * tmpx = new double [INI.Dimension()];


    	int size = 0;
		MPI_Comm_size (MPI_COMM_WORLD, &size);


    	// List of the simulations (and seeds) to be run
    	std::list<std::pair<SearchSpElement *,unsigned int> > SimulationList;

    	// List with all the available processes
    	std::list<unsigned int> AvailableProcesses;

    	// Create the list of available processes
		for (unsigned int i=1; i<(unsigned int) size; ++i){
			AvailableProcesses.push_back(i);
		}

    	// List with the species being simulated in each processor
    	std::vector<SearchSpElement *> SimulatingVector(size, (SearchSpElement *) 0);

    	// Map associating the SpElements with their associated results
    	std::map<SearchSpElement *, std::list<double> > SpElementMap;


    	// Non-blocking reception to know when a simulation is finished
    	MPI_Irecv(&result, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);

    	while(endedProcesses<size){

    		// Check if a simulation is finished
    		MPI_Test(&request, &simulation_finished, &status);

    		if(simulation_finished==1 || (AvailableProcesses.empty())){
    			if (!simulation_finished){
    				// Waiting for some process to finish
    				MPI_Wait(&request, &status);
    			}

    			// Get the species just finished and add the simulation result
    			SearchSpElement * element_finished = SimulatingVector[status.MPI_SOURCE];
    			SpElementMap[element_finished].push_back(result);

    			if (SpElementMap[element_finished].size()==INI.NumberOfSimulations()){
    				// All the simulations of this species are finished
    				double Acum = 0.0;
    				while (!SpElementMap[element_finished].empty()){
    					Acum += SpElementMap[element_finished].front();
    					SpElementMap[element_finished].pop_front();
    				}
    				element_finished->SetValue(Acum/(double)INI.NumberOfSimulations());
    				// Add this element to the output queue
    				output_queue.add(element_finished);
    				SpElementMap.erase(element_finished);
    			}

    			// Mark this process as available
    			AvailableProcesses.push_back(status.MPI_SOURCE);
    			SimulatingVector[status.MPI_SOURCE] = 0;

    			// Create a new recv request
    			MPI_Irecv(&result, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
    		} else if (not AvailableProcesses.empty()){
    			// There are available processes
    			if (not SimulationList.empty()){
    				// There are simulations ready to run
    				std::pair<SearchSpElement *, unsigned int> new_sim = SimulationList.front();
    				SimulationList.pop_front();
    				SearchSpElement * parameters = new_sim.first;
    				unsigned int seed = new_sim.second;

    				// Pack the parameters and the seed to send
    				int position=0;
    				parameters->GetX(tmpx);

    				MPI_Pack(tmpx, INI.Dimension(), MPI_DOUBLE, sendBuf, Nbytes, &position, MPI_COMM_WORLD);
    				MPI_Pack(&seed, 1, MPI_UNSIGNED, sendBuf, Nbytes,&position,MPI_COMM_WORLD);

    				// Send the simulation to one of the available processes
    				unsigned int ProcId = AvailableProcesses.front();
    				AvailableProcesses.pop_front();
    				// Assign the processor to the parameters
					SimulatingVector[ProcId] = parameters;

					if(MPI_Send(sendBuf, position, MPI_PACKED, ProcId, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    					printf("An error in MPI_Send() from master to slave\n");
    			} else if (input_queue.size()>0) {
    				// Extract a new parameter from the input queue
    				SearchSpElement * parameters = input_queue.remove();

    				// Add the parameters to the map
					if (SpElementMap.find(parameters)==SpElementMap.end()){
						SpElementMap[parameters] = std::list<double>();
					}

    				for (unsigned int i=0; i<INI.NumberOfSimulations(); ++i){
    					SimulationList.push_back(std::pair<SearchSpElement *, unsigned int>(parameters, INI.SimulationSeed()+i));
    				}
    			} else if (end_simulation){
    				// Send an ending signal to one of the available process
					unsigned int ProcId = AvailableProcesses.front();
					AvailableProcesses.pop_front();

					if(MPI_Send(sendBuf, 0, MPI_PACKED, ProcId, 1, MPI_COMM_WORLD) != MPI_SUCCESS)
						printf("An error in MPI_Send() of ending signal from master to slave\n");

					endedProcesses++;
    			} else {
    				// There is nothing to do. Sleep for 1s
    				struct timespec tim;
    				tim.tv_sec = 0;
    				tim.tv_nsec = (long)1e8;
    				nanosleep(&tim, NULL);
    			}
    		} else {
    			// There is nothing to do. Sleep for 1s
    			struct timespec tim;
				tim.tv_sec = 0;
				tim.tv_nsec = (long)1e8;
				nanosleep(&tim, NULL);
    		}
    	}

    	this->join();

        return NULL;
    };
};

#endif
