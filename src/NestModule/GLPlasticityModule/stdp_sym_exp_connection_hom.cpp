/*
 *  stdp_sym_exp_connection_hom.cpp
 *
 *  This file has been reimplemented from stdp_sym_connection_hom.h which is
 *  distributed as part of NEST.
 *
 */

#include "network.h"
#include "dictdatum.h"
#include "connector_model.h"
#include "common_synapse_properties.h"
#include "stdp_sym_exp_connection_hom.h"
#include "event.h"

namespace mynest
{
  //
  // Implementation of class STDPHomCommonProperties.
  //

  STDPSymExpHomCommonProperties::STDPSymExpHomCommonProperties() :
    CommonSynapseProperties(),
    tau_plus_(20.0),
    lambda_(0.01),
    alpha_(1.0),
    Wmax_(100.0)
  { }

  void STDPSymExpHomCommonProperties::get_status(DictionaryDatum & d) const
  {
    CommonSynapseProperties::get_status(d);

    def<nest::double_t>(d, "tau_plus", tau_plus_);
    def<nest::double_t>(d, "lambda", lambda_);
    def<nest::double_t>(d, "alpha", alpha_);
    def<nest::double_t>(d, "Wmax", Wmax_);
  }
  
  void STDPSymExpHomCommonProperties::set_status(const DictionaryDatum & d, nest::ConnectorModel &cm)
  {
    CommonSynapseProperties::set_status(d, cm);

    updateValue<nest::double_t>(d, "tau_plus", tau_plus_);
    updateValue<nest::double_t>(d, "lambda", lambda_);
    updateValue<nest::double_t>(d, "alpha", alpha_);
    updateValue<nest::double_t>(d, "Wmax", Wmax_);
  }


  //
  // Implementation of class STDPConnectionHom.
  //

  STDPSymExpConnectionHom::STDPSymExpConnectionHom() :
    Kplus_(0.0)
  { }

  STDPSymExpConnectionHom::STDPSymExpConnectionHom(const STDPSymExpConnectionHom &rhs) :
    ConnectionHetWD(rhs)
  {
    Kplus_ = rhs.Kplus_;
  }

  void STDPSymExpConnectionHom::get_status(DictionaryDatum & d) const
  {

    // base class properties, different for individual synapse
    ConnectionHetWD::get_status(d);

    // own properties, different for individual synapse
    def<nest::double_t>(d, "Kplus", Kplus_);
  }
  
  void STDPSymExpConnectionHom::set_status(const DictionaryDatum & d, nest::ConnectorModel &cm)
  {
    // base class properties
    ConnectionHetWD::set_status(d, cm);
    updateValue<nest::double_t>(d, "Kplus", Kplus_);

    // exception throwing must happen after setting own parameters
    // this exception will be caught and ignored within generic_connector_model::set_status()
    // it will not be caught if set_status is called directly to signify the specified error 
    if (d->known("tau_plus") ||
	d->known("lambda") ||
	d->known("alpha") ||
	d->known("Wmax") )
      {
	//cm.network().message(SLIInterpreter::M_ERROR, "STDPConnectionHom::set_status()", "you are trying to set common properties via an individual synapse.");
	//throw BadParameter();
	throw nest::ChangeCommonPropsByIndividual("STDPSymExpConnectionHom::set_status(): you are trying to set common properties via an individual synapse.");
      }

  }

   /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  void STDPSymExpConnectionHom::set_status(const DictionaryDatum & d, nest::index p, nest::ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);
    nest::set_property<nest::double_t>(d, "Kpluss", p, Kplus_);
  
     if (d->known("tau_pluss") ||
         d->known("lambdas") ||
         d->known("alphas") ||
         d->known("Wmaxs") )
     {
       //cm.network().message(SLIInterpreter::M_ERROR, "STDPConnectionHom::set_status()", "you are trying to set common properties via an individual synapse.");
       throw nest::ChangeCommonPropsByIndividual("STDPSymExpConnectionHom::set_status(): you are trying to set common properties via an individual synapse.");
     }

  }

  void STDPSymExpConnectionHom::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);
    initialize_property_array(d, "Kpluss");
  }
  
  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void STDPSymExpConnectionHom::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);
    append_property<nest::double_t>(d, "Kpluss", Kplus_);
  }

} // of namespace nest
