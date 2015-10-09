/*
 *  stdp_sym_exp_connection_hom.h
 *
 *  This file is based on the stdp_sym_connection_hom.h included as part of NEST.
 */

#ifndef ESTDP_CONNECTION_HOM_H
#define ESTDP_CONNECTION_HOM_H

/* BeginDocumentation
  Name: estdp_synapse_hom - Synapse type for classical exponential spike-timing dependent
   plasticity for excitatory synapses using homogeneous parameters,
   i.e. all synapses have the same parameters.

  Description:
   estdp_synapse_hom is a connector to create synapses with additive spike time
   dependent plasticity. This is a simplified version of the STDP version included with nest
   only for additive connections.
   
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

#include "connection.h"
#include "archiving_node.h"
#include <cmath>

namespace mynest
{

  /**
   * Class containing the common properties for all synapses of type STDPConnectionHom.
   */
  class ESTDPHomCommonProperties : public nest::CommonSynapseProperties
    {
    public:

      /**
       * Default constructor.
       * Sets all property values to defaults.
       */
      ESTDPHomCommonProperties():CommonSynapseProperties(),
	    tau_plus_(20.0),
	    lambda_(0.01),
	    alpha_(1.0),
	    Wmax_(100.0)
	  { };
   
      /**
       * Get all properties and put them into a dictionary.
       */
      void get_status(DictionaryDatum & d) const
      {
        CommonSynapseProperties::get_status(d);

        def<nest::double_t>(d, "tau_plus", tau_plus_);
        def<nest::double_t>(d, "lambda", lambda_);
        def<nest::double_t>(d, "alpha", alpha_);
        def<nest::double_t>(d, "Wmax", Wmax_);
      };
  
      /**
       * Set properties from the values given in dictionary.
       */
      void set_status(const DictionaryDatum & d, nest::ConnectorModel& cm)
      {
        CommonSynapseProperties::set_status(d, cm);

        updateValue<nest::double_t>(d, "tau_plus", tau_plus_);
        updateValue<nest::double_t>(d, "lambda", lambda_);
        updateValue<nest::double_t>(d, "alpha", alpha_);
        updateValue<nest::double_t>(d, "Wmax", Wmax_);
      };

      // data members common to all connections
      nest::double_t tau_plus_;
      nest::double_t lambda_;
      nest::double_t alpha_;
      nest::double_t Wmax_;
  };

  /**
   * Class representing an STDP connection with homogeneous parameters, i.e. parameters are the same for all synapses.
   */
  template < typename targetidentifierT >
  class ESTDPConnectionHom : public nest::Connection<targetidentifierT>
  {

  public:

  typedef ESTDPHomCommonProperties CommonPropertiesType;
  typedef nest::Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  ESTDPConnectionHom();
  
  /**
   * Copy constructor from a property object.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  ESTDPConnectionHom(const ESTDPConnectionHom &);

  // Explicitly declare all methods inherited from the dependent base ConnectionBase.
  // This avoids explicit name prefixes in all places these functions are used.
  // Since ConnectionBase depends on the template parameter, they are not automatically
  // found in the base class.
  using ConnectionBase::get_delay;
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;

  /**
   * Default Destructor.
   */
  virtual ~ESTDPConnectionHom() {}

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  void get_status(DictionaryDatum & d) const;
  
  /**
   * Set properties of this connection from the values given in dictionary.
   */
  void set_status(const DictionaryDatum & d, nest::ConnectorModel &cm);

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param t_lastspike Point in time of last spike sent.
   */
  void send(nest::Event& e, nest::thread t, nest::double_t t_lastspike, const ESTDPHomCommonProperties &);

  void set_weight( double_t w )
    {
      weight_ = w;
    }

  class ConnTestDummyNode : public nest::ConnTestDummyNodeBase
  {
  public:
	// Ensure proper overriding of overloaded virtual functions.
	// Return values from functions are ignored.
	using ConnTestDummyNodeBase::handles_test_event;
	nest::port handles_test_event( nest::SpikeEvent&, nest::rport )
	{
	  return nest::invalid_port_;
	}
  };

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
  void check_connection( nest::Node& s, nest::Node& t, nest::rport receptor_type, nest::double_t t_lastspike, const CommonPropertiesType& ){
    ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

    t.register_stdp_connection( t_lastspike - get_delay() );
  }
  
 private:

  nest::double_t postsynaptic_change_(nest::double_t w, nest::double_t kplus, const ESTDPHomCommonProperties &cp);
  nest::double_t presynaptic_change_(nest::double_t w, nest::double_t kminus, const ESTDPHomCommonProperties &cp);

  // data members of each connection
  nest::double_t weight_;
  nest::double_t Kplus_;

  };

template < typename targetidentifierT >
ESTDPConnectionHom< targetidentifierT >::ESTDPConnectionHom()
    : ConnectionBase()
    , weight_( 1.0 )
    , Kplus_( 0.0 )
{
}

template < typename targetidentifierT >
ESTDPConnectionHom< targetidentifierT >::ESTDPConnectionHom( const ESTDPConnectionHom& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  , Kplus_( rhs.Kplus_ )
{
}

template < typename targetidentifierT >
inline nest::double_t ESTDPConnectionHom< targetidentifierT >::postsynaptic_change_(nest::double_t w, nest::double_t kplus, const ESTDPHomCommonProperties &cp)
{
	nest::double_t norm_w = (w / cp.Wmax_) + (cp.lambda_ * kplus);
	if (norm_w > 1.0){
		return cp.Wmax_;
	}else if (norm_w < 0.0)
		return 0.0;

	return norm_w * cp.Wmax_;
}

template < typename targetidentifierT >
inline nest::double_t ESTDPConnectionHom< targetidentifierT >::presynaptic_change_(nest::double_t w, nest::double_t kminus, const ESTDPHomCommonProperties &cp)
{
	nest::double_t norm_w = (w / cp.Wmax_) - (cp.lambda_ * cp.alpha_ * kminus);
	if (norm_w > 1.0){
		return cp.Wmax_;
	}else if (norm_w < 0.0)
		return 0.0;

	return norm_w * cp.Wmax_;
}

/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
template < typename targetidentifierT >
inline void ESTDPConnectionHom< targetidentifierT >::send(nest::Event& e, nest::thread t, nest::double_t t_lastspike, const ESTDPHomCommonProperties &cp)
{
  // synapse STDP depressing/facilitation dynamics

	nest::double_t t_spike = e.get_stamp().get_ms();

  
  // t_lastspike_ = 0 initially
  
	nest::Node* target = get_target( t );
	nest::double_t dendritic_delay = nest::Time(nest::Time::step(get_delay())).get_ms();
    
  //get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque<nest::histentry>::iterator start;
  std::deque<nest::histentry>::iterator finish;
  target->get_history(t_lastspike - dendritic_delay, t_spike - dendritic_delay,
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
  weight_ = presynaptic_change_(weight_, target->get_K_value(t_spike - dendritic_delay), cp);
  e.set_receiver(*target);
  e.set_weight(weight_);
  e.set_delay(get_delay_steps());
  e.set_rport(get_rport());
  e();

  Kplus_ = Kplus_ * std::exp((t_lastspike - t_spike) /  cp.tau_plus_) + 1.0;
}

template < typename targetidentifierT >
void ESTDPConnectionHom< targetidentifierT >::get_status( DictionaryDatum& d ) const
{

  // base class properties, different for individual synapse
  ConnectionBase::get_status( d );
  def< nest::double_t >( d, nest::names::weight, weight_ );

  // own properties, different for individual synapse
  def< nest::double_t >( d, "Kplus", Kplus_ );
  def< nest::long_t >( d, nest::names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void ESTDPConnectionHom< targetidentifierT >::set_status( const DictionaryDatum& d, nest::ConnectorModel& cm )
{
  // base class properties
  ConnectionBase::set_status( d, cm );
  updateValue< nest::double_t >( d, nest::names::weight, weight_ );

  updateValue< nest::double_t >( d, "Kplus", Kplus_ );

  // exception throwing must happen after setting own parameters
  // this exception will be caught and ignored within generic_connector_model::set_status()
  // it will not be caught if set_status is called directly to signify the specified error
  // if (d->known("tau_plus") || d->known("lambda") || d->known("alpha") || d->known("Wmax") ){
  //	  throw ChangeCommonPropsByIndividual("ESTDPConnectionHom::set_status(): you are trying to set common properties via an individual synapse.");
  // }
}

} // of namespace nest

#endif // of #ifndef STDP_CONNECTION_HOM_H
