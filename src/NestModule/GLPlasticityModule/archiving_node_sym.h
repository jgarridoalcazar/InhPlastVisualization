/*
 *  Archiving_Node_Sym.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/**
 * \file Archiving_Node_Sym.h
 * Definition of Archiving_Node which is capable of
 * recording and managing a spike history.
 * \author Moritz Helias, Abigail Morrison
 * \date april 2006
 */

#ifndef ARCHIVING_NODE_SYM_H
#define ARCHIVING_NODE_SYM_H

#include "nest.h"
#include "archiving_node.h"
#include "dictdatum.h"
#include "nest_time.h"
#include "histentry_sym.h"
#include <deque>

#define DEBUG_ARCHIVER 1

namespace nest {

  namespace names
	{
    	// Neuron parameters
    	extern const Name tau_sym1;      	//!< Distance between the central and the external peaks in ms
    }
}

namespace mynest {

/**
 * \class Archiving_Node
 * a node which archives spike history for the purposes of
 * timing dependent plasticity
 */
  class Archiving_Node_Sym: public nest::Archiving_Node
{
   
 public:
  

  /**
   * \fn Archiving_Node()
   * Constructor.
   */
  Archiving_Node_Sym();

  /**
   * \fn Archiving_Node_Sym()
   * Copy Constructor.
   */
  Archiving_Node_Sym(const Archiving_Node_Sym&);


  /**
   * \fn double_t get_sym_K_value(long_t t)
   * \param central The symKminus of the central part (exp*cos^2) value at t (in ms).
   * \param external The symKminos of the external part (exp*sin^2) value at t (in ms).
   */
  void get_sym_K_value(double t, double & central, double & external);

  /**
   * \fn void get_history(long_t t1, long_t t2, std::deque<Archiver::histentry_sym>::iterator* start, std::deque<Archiver::histentry_sym>::iterator* finish)
   * return the spike times (in steps) of spikes which occurred in the range (t1,t2].
   */
  void get_sym_history(double t1, double t2,
               std::deque<histentry_sym>::iterator* start,
    		   std::deque<histentry_sym>::iterator* finish);

    /**
     * Register a new incoming STDP connection.
     *
     * t_first_read: The newly registered synapse will read the history entries with t > t_first_read.
     */
    void register_stdp_connection_sym(double t_first_read);

    void get_status(DictionaryDatum & d) const;
    void set_status(const DictionaryDatum & d);

 protected:

  /**
   * \fn void set_spiketime(Time const & t_sp)
   * record spike history
   */
  void set_spiketime(nest::Time const & t_sp);

  /**
   * \fn void clear_history()
   * clear spike history
   */
  void clear_history();


 private:

  // number of incoming connections from istdp connectors.
    // needed to determine, if every incoming connection has
    // read the spikehistory for a given point in time
    size_t n_incoming_sym_;

    // Accumulation variables
    double Kexpt1_; 		// value of exp(-abs(x)/tau1) at that time
    double Kcos2t1_;		// value of exp(-abs(x)/tau1)*cos(2*x/tau1) at that time
    double Ksin2t1_;		// value of exp(-abs(x)/tau1)*sin(2*x/tau1) at that time
    double Kexpt2_; 		// value of exp(-abs(x)/tau2) at that time
    double Ksin2t2_;		// value of exp(-abs(x)/tau2)*sin(2*x/tau2) at that time
    double Kcos2t2_;		// value of exp(-abs(x)/tau1)*cos(2*x/tau2) at that time

    // time constants for symmetric STDP
    double inv_tau_sym1_; // Distance between peaks of the kernel
    double inv_tau_sym2_; // It has to be calculated from tau1 according to tau2 = tau1/(atan(pi/2)*2*pi)

    double last_spike_sym_;

    // spiking history needed by stdp synapses
    std::deque<histentry_sym> history_sym_;

};
  
} // of namespace

#endif



