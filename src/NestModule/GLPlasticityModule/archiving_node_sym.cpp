/*
 *  Archiving_Node_Sym.cpp
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
 * \file archiving_node.cpp
 * Implementation of archiving_node to record and manage spike history
 * \author Moritz Helias, Abigail Morrison
 * \date april 2006
 */

#include "archiving_node_sym.h"
#include "dictutils.h"
#include <cmath>
#include <cstdlib>

namespace nest
{
	namespace names
  {

  	  const Name tau_sym1("tau_sym");
  }
}

namespace mynest {


  //member functions for Archiving_Node

mynest::Archiving_Node_Sym::Archiving_Node_Sym() :
		n_incoming_sym_(0),
		Kexpt1_(0.0),
		Kcos2t1_(0.0),
		Ksin2t1_(0.0),
		Kexpt2_(0.0),
		Ksin2t2_(0.0),
		Kcos2t2_(0.0),
		inv_tau_sym1_(5.e-2),
		last_spike_sym_(0.0)
		{
		inv_tau_sym2_ = inv_tau_sym1_*(std::atan(M_PI/2)*2./M_PI);
		}

mynest::Archiving_Node_Sym::Archiving_Node_Sym(const Archiving_Node_Sym& n)
:Archiving_Node(n),
     n_incoming_sym_(n.n_incoming_sym_),
     Kexpt1_(n.Kexpt1_),
     Kcos2t1_(n.Kcos2t1_),
     Ksin2t1_(n.Ksin2t1_),
     Kexpt2_(n.Kexpt2_),
     Ksin2t2_(n.Ksin2t2_),
     Kcos2t2_(n.Kcos2t2_),
	   inv_tau_sym1_(n.inv_tau_sym1_),
     inv_tau_sym2_(n.inv_tau_sym2_),
     last_spike_sym_(n.last_spike_sym_)
  {}

  void Archiving_Node_Sym::register_stdp_connection_sym(nest::double_t t_first_read)
  {
    // Mark all entries in the deque, which we will not read in future as read by this input
    // input, so that we savely increment the incoming number of
    // connections afterwards without leaving spikes in the history.
    // For details see bug #218. MH 08-04-22
	  for ( std::deque<mynest::histentry_sym>::iterator runner = history_sym_.begin();
	  runner != history_sym_.end() && runner->t_ <= t_first_read;
	  ++runner)
      (runner->access_counter_)++;

    n_incoming_sym_++;
  }

  void mynest::Archiving_Node_Sym::get_sym_K_value(nest::double_t t, nest::double_t & central, nest::double_t & external)
  {
  	central = external = 0.0;
    if (history_sym_.empty()) return;
    int i = history_sym_.size() - 1;
    while (i >= 0){
    	if (t > history_sym_[i].t_){

    		nest::double_t dt = history_sym_[i].t_ - t;

    		// Calculate central = exp(-abs(dt/tau1))*cos(dt*pi/(tau1*2))^2
    		nest::double_t dt_tau1 = dt*inv_tau_sym1_;

    		nest::double_t dt_pi_tau1 = M_PI*dt_tau1;
    		nest::double_t aux_expon_tau1 = std::exp(dt_tau1);
    		nest::double_t aux_cos_tau1 = std::cos(dt_pi_tau1);
    		nest::double_t aux_sin_tau1 = std::sin(dt_pi_tau1);

    		nest::double_t exp_t1 = history_sym_[i].Kexpt1_*aux_expon_tau1;

    		// Calculate cos_2_t1 = Cos(dt*pi/tau1)
    		nest::double_t cos_2_t1 = aux_expon_tau1*(history_sym_[i].Kcos2t1_*aux_cos_tau1 -
    				history_sym_[i].Ksin2t1_*aux_sin_tau1);

    		central = 0.5*(exp_t1 + cos_2_t1);

//      		if (central<0.0){
//      			std::cout << "Error: Central<0.0. t=" << t << " exp_t1=" << exp_t1 << " and cos_2_t1=" << cos_2_t1 << std::endl;
//      		}

    		// Calculate external = A*exp(-2*abs(dt)/tau2)*sin(dt*pi/(tau2*2))^2
    		nest::double_t dt_tau2 = dt*inv_tau_sym2_;
    		nest::double_t dt_pi_tau2 = M_PI*dt_tau2;
    		nest::double_t aux_expon_tau2 = std::exp(2*dt_tau2);
    		nest::double_t aux_cos_tau2 = std::cos(dt_pi_tau2);
    		nest::double_t aux_sin_tau2 = std::sin(dt_pi_tau2);

    		nest::double_t exp_t2 = history_sym_[i].Kexpt2_*aux_expon_tau2;

			// Calculate cos_2_t2 = Cos(dt*pi/tau2)
			nest::double_t cos_2_t2 = aux_expon_tau2*(history_sym_[i].Kcos2t2_*aux_cos_tau2 +
					history_sym_[i].Ksin2t2_*aux_sin_tau2);

			external = 0.5*(exp_t2 - cos_2_t2);

			//std::cout << "Post-pre traces. Central: " << central << " Expt1: " << exp_t1 << " Cos2t1: " << cos_2_t1 << " External: " << external << " Expt2: " << exp_t2 << " Cos2t2: " << cos_2_t2 << std::endl;


//  			if (external<0.0){
//				std::cout << "Error: External<0.0. t=" << t << " exp_t2=" << exp_t2 << " and cos_2_t2=" << cos_2_t2 << std::endl;
//			}

//  			std::cout << "Node trace: central: " << central << " external: " << external << std::endl;
    		return;
    	}
    	i--;
    }
    return;
  }

  void mynest::Archiving_Node_Sym::get_sym_history(nest::double_t t1, nest::double_t t2,
  				   std::deque<histentry_sym>::iterator* start,
  				   std::deque<histentry_sym>::iterator* finish)
  {

	  *finish = history_sym_.end();
      if (history_sym_.empty())
        {
  	*start = *finish;
  	return;
        }
      else
        {
  	std::deque<histentry_sym>::iterator runner = history_sym_.begin();
  	while ((runner != history_sym_.end()) && (runner->t_ <= t1)) ++runner;
  	*start = runner;
  	while ((runner != history_sym_.end()) && (runner->t_ <= t2))
  	  {
  	    (runner->access_counter_)++;
  	    ++runner;
  	  }
  	*finish = runner;
        }
  }

  void mynest::Archiving_Node_Sym::set_spiketime(nest::Time const & t_sp)
  {
    Archiving_Node::set_spiketime(t_sp);

    if (n_incoming_sym_){
      // prune all spikes from history which are no longer needed
      // except the penultimate one. we might still need it.
      while (history_sym_.size() > 1){
        if (history_sym_.front().access_counter_ >= n_incoming_sym_)
          history_sym_.pop_front();
        else
          break;
      }      

      //std::cout << "Updating postsynaptic neuron from time " << last_spike_sym_ << " to " << t_sp.get_ms() << std::endl;
      nest::double_t dt = last_spike_sym_ - t_sp.get_ms();

      // Calculate central = exp(-abs(dt/tau1))*cos(dt*pi/(tau1*2))^2
      nest::double_t dt_tau1 = dt*inv_tau_sym1_;

      nest::double_t dt_pi_tau1 = M_PI*dt_tau1;
      nest::double_t aux_expon_tau1 = std::exp(dt_tau1);
      nest::double_t aux_cos_tau1 = std::cos(dt_pi_tau1);
      nest::double_t aux_sin_tau1 = std::sin(dt_pi_tau1);

      Kexpt1_ = Kexpt1_*aux_expon_tau1 + 1.0;

      // Calculate cos_2_t1 = Cos(dt*pi/tau1)
      Kcos2t1_ = aux_expon_tau1*(Kcos2t1_*aux_cos_tau1 - Ksin2t1_*aux_sin_tau1) + 1.0;
      Ksin2t1_ = aux_expon_tau1*(Ksin2t1_*aux_cos_tau1 + Kcos2t1_*aux_sin_tau1);

      // Calculate external = A*exp(-2*abs(dt)/tau2)*sin(dt*pi/(tau2*2))^2
      nest::double_t dt_tau2 = dt*inv_tau_sym2_;

      // std::cout << "Postsynaptic spike at time " << t_sp.get_ms() << " dt: " << dt << " - Kexpt1: " << Kexpt1_ << " - TauSym1: " << 1./inv_tau_sym1_ << std::endl;

      nest::double_t dt_pi_tau2 = M_PI*dt_tau2;
      nest::double_t aux_expon_tau2 = std::exp(2*dt_tau2);
      nest::double_t aux_cos_tau2 = std::cos(dt_pi_tau2);
      nest::double_t aux_sin_tau2 = std::sin(dt_pi_tau2);

      Kexpt2_ = Kexpt2_*aux_expon_tau2 + 1.0;

      // Calculate cos_2_t2 = Cos(dt*pi/tau2)
      Kcos2t2_ = aux_expon_tau2*(Kcos2t2_*aux_cos_tau2 + Ksin2t2_*aux_sin_tau2) + 1.0;
      Ksin2t2_ = aux_expon_tau2*(Ksin2t2_*aux_cos_tau2 + Kcos2t2_*aux_sin_tau2);

      last_spike_sym_ = t_sp.get_ms();

      //std::cout << "Updated postsynaptic trace. Expt1: " << Kexpt1_ << " Cos2t1: " << Kcos2t1_ << "Sin2t1: " << Ksin2t1_ << " Expt2: " << Kexpt2_ << " Cos2t2: " << Kcos2t2_ << " Sin2t2: " << Ksin2t2_ << std::endl;

      //std::cout << "Post-spike received at " << last_spike_sym_ << "- Kexpt1 - " << Kexpt1_ << " Kcos2t1 - " << Kcos2t1_ << " Ksin2t1 - " << Ksin2t1_
      // << " Kexpt2 - " << Kexpt2_ << " Ksin2t2 - " << Ksin2t2_ << " Kcos2t2 - " << Kcos2t2_ << std::endl;

      history_sym_.push_back( histentry_sym( last_spike_sym_, Kexpt1_, Kcos2t1_, Ksin2t1_, Kexpt2_, Ksin2t2_, Kcos2t2_, 0) );
    } else {
      last_spike_sym_ = t_sp.get_ms();
    }
  }


  void mynest::Archiving_Node_Sym::get_status(DictionaryDatum & d) const
  {
	  Archiving_Node::get_status(d);
    def<nest::double_t>(d, nest::names::t_spike, get_spiketime_ms());
    def<nest::double_t>(d, nest::names::tau_sym1, 1.0/inv_tau_sym1_);
  #ifdef DEBUG_ARCHIVER
    def<int>(d, nest::names::archiver_length, history_sym_.size());
  #endif
  }

  void mynest::Archiving_Node_Sym::set_status(const DictionaryDatum & d)
  {
	  Archiving_Node::set_status(d);
    // We need to preserve values in case invalid values are set
	  nest::double_t new_tau_sym1 = 1.0/inv_tau_sym1_;
    updateValue<nest::double_t>(d, nest::names::tau_sym1, new_tau_sym1);

    if ( new_tau_sym1 <= 0)
      throw nest::BadProperty("All time constants must be strictly positive.");

    inv_tau_sym1_ = 1.0/new_tau_sym1;
    inv_tau_sym2_ = inv_tau_sym1_*(std::atan(M_PI/2)*2./M_PI);

    // check, if to clear spike history and K_minus
    bool clear = false;
    updateValue<bool>(d, nest::names::clear, clear);
    if ( clear )
    	clear_history();
  }

  void mynest::Archiving_Node_Sym::clear_history()
  {
  	Archiving_Node::clear_history();

  	last_spike_sym_ = 0.0;
    Kexpt1_ = 0.0;
    Kcos2t1_ = 0.0;
    Ksin2t1_ = 0.0;
   	Kexpt2_ = 0.0;
   	Kcos2t2_ = 0.0;
    Ksin2t2_ = 0.0;
    history_sym_.clear();
  }

} // of namespace nest
