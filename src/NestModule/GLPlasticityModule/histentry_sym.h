/*
 *  histentry_sym.h
 */

/**
 * \file histentry_sym.h
 * Part of definition of Archiving_Node_sym which is capable of
 * recording and managing a spike history. 
 * \author Jesus Garrido
 * \note modified from histentry.h
 * \date september 2014
 */

#ifndef HISTENTRYSYM_H
#define HISTENTRYSYM_H

#include "nest.h"

namespace mynest {

// entry in the spiking history
  class histentry_sym 
  {
    public:
      histentry_sym(nest::double_t t, nest::double_t Kexpt1, 
        nest::double_t Kcos2t1, nest::double_t Ksin2t1,
			  nest::double_t Kexpt2, nest::double_t Ksin2t2, nest::double_t Kcos2t2,
			  size_t access_counter);

      nest::double_t t_;              // point in time when spike occurred (in ms)
      size_t access_counter_;   // how often this entry was accessed (to enable removal,
  
      nest::double_t Kexpt1_; 		// value of exp(-abs(x)/tau1) at that time
      nest::double_t Kcos2t1_;		// value of cos(2*x/tau1) at that time
      nest::double_t Ksin2t1_;		// value of sin(2*x/tau1) at that time
      nest::double_t Kexpt2_; 		// value of exp(-abs(x)/tau2) at that time
      nest::double_t Ksin2t2_;		// value of sin(2*x/tau2) at that time
      nest::double_t Kcos2t2_;		// value of cos(2*x/tau2) at that time

  };

}

#endif
