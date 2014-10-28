/*
 *  histentry_sym.cpp
 */

/**
 * \file histentry_sym.cpp
 * Implementation of archiving_node to record and manage spike history
 * \author Jesus Garrido
 * \date september 2014
 * \note reimplemented from histentry.cpp
 */

#include "histentry_sym.h"

  // member functions of histentry_sym

  mynest::histentry_sym::histentry_sym(nest::double_t t, nest::double_t Kminus, nest::double_t triplet_Kminus,
		  nest::double_t Kexpt1, nest::double_t Kcos2t1, nest::double_t Ksin2t1, nest::double_t Kexpt2,
		  nest::double_t Ksin2t2, nest::double_t Kcos2t2, size_t access_counter) :
    nest::histentry(t, Kminus, triplet_Kminus, access_counter),
	Kexpt1_(Kexpt1),
	Kcos2t1_(Kcos2t1),
	Ksin2t1_(Ksin2t1),
	Kexpt2_(Kexpt2),
	Ksin2t2_(Ksin2t2),
	Kcos2t2_(Kcos2t2)
  { } 

