/*
 *  stdp_sym_connection_hom.cpp
 *
 *  This file has been reimplemented from stdp_sym_connection_hom.h which is
 *  distributed as part of NEST.
 *
 */

#include "network.h"
#include "dictdatum.h"
#include "connector_model.h"
#include "common_synapse_properties.h"
#include "stdp_sym_connection_hom.h"
#include "event.h"

namespace mynest
{
  //
  // Implementation of class STDPSymHomCommonProperties.
  //

  STDPSymHomCommonProperties::STDPSymHomCommonProperties() :
    CommonSynapseProperties(),
    inv_tau_sym1_(50.e-3),
    lambda_(0.01),
    alpha_(0.25),
    Wmax_(100.0)
  {
    inv_tau_sym2_ = inv_tau_sym1_*(std::atan(M_PI/2)*2./M_PI);

    sym_A_ = 1./(std::exp(-std::atan(M_PI/2.)*4./M_PI)*std::pow(std::sin(std::atan(M_PI/2.)),2)); // Normalizing constant of the external part
  }

  void STDPSymHomCommonProperties::get_status(DictionaryDatum & d) const
  {
    CommonSynapseProperties::get_status(d);

    def<nest::double_t>(d, "tau_sym", 1./inv_tau_sym1_);
    def<nest::double_t>(d, "lambda", lambda_);
    def<nest::double_t>(d, "alpha", alpha_);
    def<nest::double_t>(d, "Wmax", Wmax_);
  }
  
  void STDPSymHomCommonProperties::set_status(const DictionaryDatum & d, nest::ConnectorModel &cm)
  {
    CommonSynapseProperties::set_status(d, cm);

    nest::double_t old_tau_sym = 1.0/inv_tau_sym1_;

    updateValue<nest::double_t>(d, "tau_sym", old_tau_sym);

    if ( old_tau_sym <= 0)
    	throw nest::BadProperty("All time constants must be strictly positive.");

    inv_tau_sym1_ = 1.0/old_tau_sym;
    inv_tau_sym2_ = inv_tau_sym1_*(std::atan(M_PI/2)*2./M_PI);

    updateValue<nest::double_t>(d, "lambda", lambda_);
    updateValue<nest::double_t>(d, "alpha", alpha_);
    updateValue<nest::double_t>(d, "Wmax", Wmax_);
  }


  //
  // Implementation of class STDPSymConnectionHom.
  //

  STDPSymConnectionHom::STDPSymConnectionHom() :
      Kexpt1_(0.0),
      Kcos2t1_(0.0),
      Ksin2t1_(0.0),
      Kexpt2_(0.0),
      Ksin2t2_(0.0),
      Kcos2t2_(0.0)
  { }

  STDPSymConnectionHom::STDPSymConnectionHom(const STDPSymConnectionHom &rhs) :
    ConnectionHetWD(rhs)
  {
	  Kexpt1_ = rhs.Kexpt1_;
	  Kcos2t1_ = Kcos2t1_;
	  Ksin2t1_ = Ksin2t1_;
	  Kexpt2_ = Kexpt2_;
	  Ksin2t2_ = Ksin2t2_;
	  Kcos2t2_ = Kcos2t2_;
  }

  void STDPSymConnectionHom::get_status(DictionaryDatum & d) const
  {

    // base class properties, different for individual synapse
    ConnectionHetWD::get_status(d);

    // own properties, different for individual synapse
    def<nest::double_t>(d, "Kexpt1", Kexpt1_);
    def<nest::double_t>(d, "Kcos2t1", Kcos2t1_);
    def<nest::double_t>(d, "Ksin2t1", Ksin2t1_);
    def<nest::double_t>(d, "Kexpt2", Kexpt2_);
    def<nest::double_t>(d, "Kcos2t2", Kcos2t2_);
    def<nest::double_t>(d, "Ksin2t2", Ksin2t2_);
  }
  
  void STDPSymConnectionHom::set_status(const DictionaryDatum & d, nest::ConnectorModel &cm)
  {
    // base class properties
    ConnectionHetWD::set_status(d, cm);
    updateValue<nest::double_t>(d, "Kexpt1", Kexpt1_);
    updateValue<nest::double_t>(d, "Kcos2t1", Kcos2t1_);
    updateValue<nest::double_t>(d, "Ksin2t1", Ksin2t1_);
    updateValue<nest::double_t>(d, "Kexpt2", Kexpt2_);
    updateValue<nest::double_t>(d, "Kcos2t2", Kcos2t2_);
    updateValue<nest::double_t>(d, "Ksin2t2", Ksin2t2_);

    // exception throwing must happen after setting own parameters
    // this exception will be caught and ignored within generic_connector_model::set_status()
    // it will not be caught if set_status is called directly to signify the specified error 
    if (d->known("tau_sym") ||
	d->known("lambda") ||
	d->known("alpha") ||
	d->known("Wmax") )
      {
	//cm.network().message(SLIInterpreter::M_ERROR, "STDPConnectionHom::set_status()", "you are trying to set common properties via an individual synapse.");
	//throw BadParameter();
	throw nest::ChangeCommonPropsByIndividual("STDPSymConnectionHom::set_status(): you are trying to set common properties via an individual synapse.");
      }

  }

   /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  void STDPSymConnectionHom::set_status(const DictionaryDatum & d, nest::index p, nest::ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);
    nest::set_property<nest::double_t>(d, "Kexpt1", p, Kexpt1_);
    nest::set_property<nest::double_t>(d, "Kcos2t1", p, Kcos2t1_);
    nest::set_property<nest::double_t>(d, "Ksin2t1", p, Ksin2t1_);
    nest::set_property<nest::double_t>(d, "Kexpt2", p, Kexpt2_);
    nest::set_property<nest::double_t>(d, "Kcos2t2", p, Kcos2t2_);
    nest::set_property<nest::double_t>(d, "Ksin2t2", p, Ksin2t2_);
  
     if (d->known("tau_syms") ||
         d->known("lambdas") ||
         d->known("alphas") ||
         d->known("Wmaxs") )
     {
       //cm.network().message(SLIInterpreter::M_ERROR, "STDPConnectionHom::set_status()", "you are trying to set common properties via an individual synapse.");
       throw nest::ChangeCommonPropsByIndividual("STDPSymConnectionHom::set_status(): you are trying to set common properties via an individual synapse.");
     }

  }

  void STDPSymConnectionHom::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);
    initialize_property_array(d, "Kexpt1s");
    initialize_property_array(d, "Kcos2t1s");
    initialize_property_array(d, "Ksin2t1s");
    initialize_property_array(d, "Kexpt2s");
    initialize_property_array(d, "Kcos2t2s");
    initialize_property_array(d, "Ksin2t2s");
  }
  
  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void STDPSymConnectionHom::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);
    append_property(d, "Kexpt1s", Kexpt1_);
    append_property(d, "Kcos2t1s", Kcos2t1_);
    append_property(d, "Ksin2t1s", Ksin2t1_);
    append_property(d, "Kexpt2s", Kexpt2_);
    append_property(d, "Kcos2t2s", Kcos2t2_);
    append_property(d, "Ksin2t2s", Ksin2t2_);
  }

} // of namespace nest
