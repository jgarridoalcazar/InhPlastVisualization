/*
 *  istdp_connection_hom.h
 *
 *  This file is based on the stdp_sym_connection_hom.h included as part of NEST.
 */

#ifndef ISTDP_CONNECTION_HOM_H
#define ISTDP_CONNECTION_HOM_H

/* BeginDocumentation
  Name: istdp_synapse_hom - Synapse type for inhibitory spike-timing dependent
   plasticity with exponential LTP and fixed LTD for each presynaptic spike.
   It uses homogeneous parameters, i.e. all synapses have the same parameters.

  Description:
   istdp_synapse_hom is a connector to create synapses with additive inhibitory spike time
   dependent plasticity.
   
  Parameters:
   tau_plus   double - Time constant of STDP window, potentiation in ms
                       (tau_minus_triplet has to be defined in post-synaptic neuron too with the same value)
   lambda     double - Step size of the potentiation at its peak.
   alpha      double - Ratio between depression and potentiation (scales depressing increments as alpha*lambda)
   Wmax       double - Maximum allowed weight


  Transmits: SpikeEvent
   
  References:
   [1] Vogels, T.P., Sprekeler, H., Zenke, F., Clopath, C., and Gerstner, W. (2011).
   	   Inhibitory plasticity balances exitation and inhibition in sensory pathways
   	   and memory networks, science 334, 1569--1573

  FirstVersion: April 2015
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
  class ISTDPHomCommonProperties : public nest::CommonSynapseProperties
    {
    public:

      /**
       * Default constructor.
       * Sets all property values to defaults.
       */
      ISTDPHomCommonProperties(): CommonSynapseProperties(),
    	tau_plus_(20.0),
    	lambda_(0.01),
    	alpha_(1.0),
    	Wmax_(100.0)
      { };
   
      /**
       * Get all properties and put them into a dictionary.
       */
      void get_status(DictionaryDatum & d) const{
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
  class ISTDPConnectionHom : public nest::Connection<targetidentifierT>
  {

  public:

  typedef ISTDPHomCommonProperties CommonPropertiesType;
  typedef nest::Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  ISTDPConnectionHom();
  
  /**
   * Copy constructor from a property object.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  ISTDPConnectionHom(const ISTDPConnectionHom &);

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
  virtual ~ISTDPConnectionHom() {}

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
  void send(nest::Event& e, nest::thread t, nest::double_t t_lastspike, const ISTDPHomCommonProperties &);

  void set_weight( double_t w )
  {
    weight_ = w;
  }

  class ConnTestDummyNode : public nest::ConnTestDummyNodeBase {
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

  nest::double_t postsynaptic_change_(nest::double_t w, nest::double_t kplus, const ISTDPHomCommonProperties &cp);
  nest::double_t presynaptic_change_(nest::double_t w, nest::double_t kminus, const ISTDPHomCommonProperties &cp);

  // data members of each connection
  nest::double_t weight_;
  nest::double_t Kplus_;

};

template < typename targetidentifierT >
ISTDPConnectionHom< targetidentifierT >::ISTDPConnectionHom()
    : ConnectionBase()
    , weight_( 1.0 )
    , Kplus_( 0.0 )
{
}

template < typename targetidentifierT >
ISTDPConnectionHom< targetidentifierT >::ISTDPConnectionHom( const ISTDPConnectionHom& rhs )
    : ConnectionBase( rhs )
    , weight_( rhs.weight_ )
    , Kplus_( rhs.Kplus_ )
{
}

template < typename targetidentifierT >
inline nest::double_t ISTDPConnectionHom< targetidentifierT >::postsynaptic_change_(nest::double_t w, nest::double_t kplus, const ISTDPHomCommonProperties &cp)
{
	nest::double_t norm_w = (w / cp.Wmax_) + (cp.lambda_ * kplus);
	if (norm_w > 1.0){
		return cp.Wmax_;
	}else if (norm_w < 0.0)
		return 0.0;

	return norm_w * cp.Wmax_;
}

template < typename targetidentifierT >
inline nest::double_t ISTDPConnectionHom< targetidentifierT >::presynaptic_change_(nest::double_t w, nest::double_t kminus, const ISTDPHomCommonProperties &cp)
{
	nest::double_t norm_w = (w / cp.Wmax_) + cp.lambda_ * (kminus - cp.alpha_);
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
inline void ISTDPConnectionHom< targetidentifierT >::send(nest::Event& e, nest::thread t, nest::double_t t_lastspike, const ISTDPHomCommonProperties &cp)
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
  double_t k_value, k_triplet;
  target->get_K_values(t_spike - dendritic_delay, k_value, k_triplet);
  weight_ = presynaptic_change_(weight_, k_triplet, cp);
  e.set_receiver(*target);
  e.set_weight(weight_);
  e.set_delay(get_delay_steps());
  e.set_rport(get_rport());
  e();

  Kplus_ = Kplus_ * std::exp((t_lastspike - t_spike) /  cp.tau_plus_) + 1.0;
  }

template < typename targetidentifierT >
void ISTDPConnectionHom< targetidentifierT >::get_status( DictionaryDatum& d ) const
{

  // base class properties, different for individual synapse
  ConnectionBase::get_status( d );
  def< nest::double_t >( d, nest::names::weight, weight_ );

  // own properties, different for individual synapse
  def< nest::double_t >( d, "Kplus", Kplus_ );
  def< nest::long_t >( d, nest::names::size_of, sizeof( *this ) );
}

} // of namespace nest

#endif // of #ifndef STDP_CONNECTION_HOM_H
