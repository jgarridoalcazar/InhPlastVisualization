/*
 *  stdp_sym_exp_connection_hom.h
 *
 *  This file is based on the stdp_sym_connection_hom.h included as part of NEST.
 */

#ifndef STDP_SYM_EXP_CONNECTION_HOM_H
#define STDP_SYM_EXP_CONNECTION_HOM_H

/* BeginDocumentation
  Name: stdp_sym_exp_synapse_hom - Synapse type for spike-timing dependent
   plasticity with y-axis symmetric kernel function using homogeneous parameters,
   i.e. all synapses have the same parameters.

  Description:
   stdp_sym_exp_synapse_hom is a connector to create synapses with additive spike time
   dependent plasticity. It produces potentiation in when presynaptic and postsynaptic activity
   is simultaneous (or close) and depression every time a presynaptic spike happens.
   
  Parameters:
   tau_plus   double - Time constant of STDP window, potentiation in ms
                       (tau_minus has to be defined in post-synaptic neuron too)
   lambda     double - Step size of the potentiation at its peak.
   alpha      double - Ratio between depression and potentiation (scales depressing increments as alpha*lambda)
   Wmax       double - Maximum allowed weight


  Transmits: SpikeEvent
   
  References:
   [2] Rubin, J., Lee, D. and Sompolinsky, H. (2001). Equilibrium
       properties of temporally asymmetric Hebbian plasticity, PRL
       86,364-367

   [3] Song, S., Miller, K. D. and Abbott, L. F. (2000). Competitive
       Hebbian learning through spike-timing-dependent synaptic
       plasticity,Nature Neuroscience 3:9,919--926

  FirstVersion: January 2015
  Author: Jesus Garrido (jesusgarrido@ugr.es)
*/

#include "connection_het_wd.h"
#include "archiving_node.h"
#include <cmath>

namespace mynest
{

  /**
   * Class containing the common properties for all synapses of type STDPConnectionHom.
   */
  class STDPSymExpHomCommonProperties : public nest::CommonSynapseProperties
    {
      friend class STDPSymExpConnectionHom;

    public:

      /**
       * Default constructor.
       * Sets all property values to defaults.
       */
      STDPSymExpHomCommonProperties();
   
      /**
       * Get all properties and put them into a dictionary.
       */
      void get_status(DictionaryDatum & d) const;
  
      /**
       * Set properties from the values given in dictionary.
       */
      void set_status(const DictionaryDatum & d, nest::ConnectorModel& cm);

      // overloaded for all supported event types
      void check_event(nest::SpikeEvent&) {}
 
    private:

      // data members common to all connections
      nest::double_t tau_plus_;
      nest::double_t lambda_;
      nest::double_t alpha_;
      nest::double_t Wmax_;
    };



  /**
   * Class representing an STDP connection with homogeneous parameters, i.e. parameters are the same for all synapses.
   */
  class STDPSymExpConnectionHom : public nest::ConnectionHetWD
  {

  public:
  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  STDPSymExpConnectionHom();
  
  /**
   * Copy constructor from a property object.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  STDPSymExpConnectionHom(const STDPSymExpConnectionHom &);

  /**
   * Default Destructor.
   */
  virtual ~STDPSymExpConnectionHom() {}

  /*
   * This function calls check_connection on the sender and checks if the receiver
   * accepts the event type and receptor type requested by the sender.
   * Node::check_connection() will either confirm the receiver port by returning
   * true or false if the connection should be ignored.
   * We have to override the base class' implementation, since for STDP
   * connections we have to call register_stdp_connection on the target neuron
   * to inform the Archiver to collect spikes for this connection.
   *
   * \param s The source node
   * \param r The target node
   * \param receptor_type The ID of the requested receptor type
   * \param t_lastspike last spike produced by presynaptic neuron (in ms)
   */
  void check_connection(nest::Node & s, nest::Node & r, nest::rport receptor_type, nest::double_t t_lastspike);

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  void get_status(DictionaryDatum & d) const;
  
  /**
   * Set properties of this connection from the values given in dictionary.
   */
  void set_status(const DictionaryDatum & d, nest::ConnectorModel &cm);

  /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  void set_status(const DictionaryDatum & d, nest::index p, nest::ConnectorModel &cm);

  /**
   * Create new empty arrays for the properties of this connection in the given
   * dictionary. It is assumed that they are not existing before.
   */
  void initialize_property_arrays(DictionaryDatum & d) const;

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void append_properties(DictionaryDatum & d) const;

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param t_lastspike Point in time of last spike sent.
   */
  void send(nest::Event& e, nest::double_t t_lastspike, const STDPSymExpHomCommonProperties &);

  // overloaded for all supported event types
  using Connection::check_event;
  void check_event(nest::SpikeEvent&) {}
  
 private:

  nest::double_t postsynaptic_change_(nest::double_t w, nest::double_t kplus, const STDPSymExpHomCommonProperties &cp);
  nest::double_t presynaptic_change_(nest::double_t w, nest::double_t kminus, const STDPSymExpHomCommonProperties &cp);

  // data members of each connection
  nest::double_t Kplus_;

  };


inline
nest::double_t STDPSymExpConnectionHom::postsynaptic_change_(nest::double_t w, nest::double_t kplus, const STDPSymExpHomCommonProperties &cp)
{
	nest::double_t norm_w = (w / cp.Wmax_) + (cp.lambda_ * kplus);
	if (norm_w > 1.0){
		return cp.Wmax_;
	}else if (norm_w < 0.0)
		return 0.0;

	return norm_w * cp.Wmax_;
}

inline
nest::double_t STDPSymExpConnectionHom::presynaptic_change_(nest::double_t w, nest::double_t kminus, const STDPSymExpHomCommonProperties &cp)
{
	nest::double_t norm_w = (w / cp.Wmax_) + (cp.lambda_ * (kminus - cp.alpha_));
	if (norm_w > 1.0){
		return cp.Wmax_;
	}else if (norm_w < 0.0)
		return 0.0;

	return norm_w * cp.Wmax_;
}


inline 
  void STDPSymExpConnectionHom::check_connection(nest::Node & s, nest::Node & r, nest::rport receptor_type, nest::double_t t_lastspike)
{
	nest::ConnectionHetWD::check_connection(s, r, receptor_type, t_lastspike);
  r.register_stdp_connection(t_lastspike - nest::Time(nest::Time::step(delay_)).get_ms());
}

/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
inline
void STDPSymExpConnectionHom::send(nest::Event& e, nest::double_t t_lastspike, const STDPSymExpHomCommonProperties &cp)
{
  // synapse STDP depressing/facilitation dynamics

	nest::double_t t_spike = e.get_stamp().get_ms();

  
  // t_lastspike_ = 0 initially
  

	nest::double_t dendritic_delay = nest::Time(nest::Time::step(delay_)).get_ms();
    
  //get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque<nest::histentry>::iterator start;
  std::deque<nest::histentry>::iterator finish;
  target_->get_history(t_lastspike - dendritic_delay, t_spike - dendritic_delay,
			       &start, &finish);
  //facilitation due to post-synaptic spikes since last pre-synaptic spike
  while (start != finish)
  {
    nest::double_t minus_dt = t_lastspike - (start->t_ + dendritic_delay);
    ++start;
    if (minus_dt == 0)
      continue;
    weight_ = postsynaptic_change_(weight_, Kplus_ * std::exp(minus_dt / cp.tau_plus_), cp);
  }


  //depression due to new pre-synaptic spike
  weight_ = presynaptic_change_(weight_, target_->get_K_value(t_spike - dendritic_delay), cp);
  e.set_receiver(*target_);
  e.set_weight(weight_);
  e.set_delay(delay_);
  e.set_rport(rport_);
  e();

  Kplus_ = Kplus_ * std::exp((t_lastspike - t_spike) /  cp.tau_plus_) + 1.0;
  }

} // of namespace nest

#endif // of #ifndef STDP_CONNECTION_HOM_H
