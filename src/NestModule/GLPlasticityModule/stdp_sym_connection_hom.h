/*
 *  stdp_sym_connection_hom.h
 *
 *  This file is based on the stdp_sym_connection_hom.h included as part of NEST.
 */

#ifndef STDP_SYM_CONNECTION_HOM_H
#define STDP_SYM_CONNECTION_HOM_H

/* BeginDocumentation
  Name: stdp_sym_synapse_hom - Synapse type for spike-timing dependent
   plasticity with y-axis symmetric kernel function using homogeneous parameters,
   i.e. all synapses have the same parameters.

  Description:
   stdp_sym_synapse_hom is a connector to create synapses with additive spike time
   dependent plasticity. It produces potentiation in when presynaptic and postsynaptic activity
   is simultaneous (or close) and depression when the happen further.
   
  Parameters:
   tau_sym   double - Distance from the central potentiation peak to the external depression peaks in ms
                       (tau_sym has to be defined in post-synaptic neuron too)
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

  FirstVersion: October 2014
  Author: Jesus Garrido (jesusgarrido@ugr.es)
*/

#include "connection.h"
#include "archiving_node_sym.h"
#include <cmath>

namespace mynest
{

  /**
   * Class containing the common properties for all synapses of type STDPConnectionHom.
   */
  class STDPSymHomCommonProperties : public nest::CommonSynapseProperties
    {
    public:

      /**
       * Default constructor.
       * Sets all property values to defaults.
       */
      STDPSymHomCommonProperties():
	    CommonSynapseProperties(),
	    inv_tau_sym1_(50.e-3),
	    lambda_(0.01),
	    alpha_(0.25),
	    Wmax_(100.0)
	  {
	    inv_tau_sym2_ = inv_tau_sym1_*(std::atan(M_PI/2)*2./M_PI);

	    sym_A_ = 1./(std::exp(-std::atan(M_PI/2.)*4./M_PI)*std::pow(std::sin(std::atan(M_PI/2.)),2)); // Normalizing constant of the external part
	  };
   
      /**
       * Get all properties and put them into a dictionary.
       */
      void get_status(DictionaryDatum & d) const
      {
    	CommonSynapseProperties::get_status(d);

    	def<double>(d, "tau_sym", 1./inv_tau_sym1_);
    	def<double>(d, "lambda", lambda_);
    	def<double>(d, "alpha", alpha_);
    	def<double>(d, "Wmax", Wmax_);
      }
  
      /**
       * Set properties from the values given in dictionary.
       */
      void set_status(const DictionaryDatum & d, nest::ConnectorModel& cm){
    	CommonSynapseProperties::set_status(d, cm);

    	double old_tau_sym = 1.0/inv_tau_sym1_;

    	updateValue<double>(d, "tau_sym", old_tau_sym);

    	if ( old_tau_sym <= 0)
    	  throw nest::BadProperty("All time constants must be strictly positive.");

    	inv_tau_sym1_ = 1.0/old_tau_sym;
    	inv_tau_sym2_ = inv_tau_sym1_*(std::atan(M_PI/2)*2./M_PI);

    	updateValue<double>(d, "lambda", lambda_);
    	updateValue<double>(d, "alpha", alpha_);
    	updateValue<double>(d, "Wmax", Wmax_);
      };

      // data members common to all connections
      double inv_tau_sym1_;
      double lambda_;
      double alpha_;
      double Wmax_;

      // Auxiliar variables common to all connections
      double inv_tau_sym2_;
      double sym_A_;
    };

  /**
   * Class representing an STDP connection with homogeneous parameters, i.e. parameters are the same for all synapses.
   */
  template < typename targetidentifierT >
  class STDPSymConnectionHom : public nest::Connection<targetidentifierT>
  {

  public:
  typedef STDPSymHomCommonProperties CommonPropertiesType;
  typedef nest::Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  STDPSymConnectionHom();
  
  /**
   * Copy constructor from a property object.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  STDPSymConnectionHom(const STDPSymConnectionHom &);

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
  virtual ~STDPSymConnectionHom() {}

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
  void send(nest::Event& e, nest::thread t, double t_lastspike, const CommonPropertiesType &);

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
  void check_connection( nest::Node& s, nest::Node& t, nest::rport receptor_type, double t_lastspike, const CommonPropertiesType& ){
	ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

    ((mynest::Archiving_Node_Sym *)(&t))->register_stdp_connection_sym( t_lastspike - get_delay() );
  }
  
 private:

  double apply_weight_change_(double w, double trace, const STDPSymHomCommonProperties &cp);

  // data members of each connection
  // Accumulation variables
  double Kexpt1_; 		// value of exp(-abs(x)/tau1) at that time
  double Kcos2t1_;		// value of exp(-abs(x)/tau1)*cos(2*x/tau1) at that time
  double Ksin2t1_;		// value of exp(-abs(x)/tau1)*sin(2*x/tau1) at that time
  double Kexpt2_; 		// value of exp(-abs(x)/tau1) at that time
  double Ksin2t2_;		// value of exp(-abs(x)/tau2)*sin(2*x/tau2) at that time
  double Kcos2t2_;		// value of exp(-abs(x)/tau2)*cos(2*x/tau2) at that time

  double weight_;
  };

template < typename targetidentifierT >
STDPSymConnectionHom< targetidentifierT >::STDPSymConnectionHom()
    : ConnectionBase(),
    weight_( 1.0 ),
	Kexpt1_(0.0),
	Kcos2t1_(0.0),
	Ksin2t1_(0.0),
	Kexpt2_(0.0),
	Ksin2t2_(0.0),
	Kcos2t2_(0.0)
{
}

template < typename targetidentifierT >
STDPSymConnectionHom< targetidentifierT >::STDPSymConnectionHom( const STDPSymConnectionHom& rhs )
    : ConnectionBase( rhs ),
	  weight_( rhs.weight_ ),
	  Kexpt1_ (rhs.Kexpt1_),
	  Kcos2t1_ (rhs.Kcos2t1_),
	  Ksin2t1_(rhs.Ksin2t1_),
      Kexpt2_(rhs.Kexpt2_),
	  Ksin2t2_(rhs.Ksin2t2_),
	  Kcos2t2_(rhs.Kcos2t2_)
{
}


template < typename targetidentifierT >
inline double STDPSymConnectionHom< targetidentifierT >::apply_weight_change_(double w, double trace, const CommonPropertiesType & cp)
{
  double norm_w = (w / cp.Wmax_) + (cp.lambda_ * trace);
//  std::cout << "Old weight: " << w << " Norm. Weight: " << norm_w << " Kexpt1: " << Kexpt1_ << " Kcos2t1: " << Kcos2t1_ << " Ksin2t1: " << Ksin2t1_ << " Kexpt2: " << Kexpt2_ << " Ksin2t2: " << Ksin2t2_ << " Kcos2t2: " << Kcos2t2_ << " Trace: " << trace << " Lambda: " << cp.lambda_ << " Wmax: " << cp.Wmax_ << this << std::endl;
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
inline void STDPSymConnectionHom< targetidentifierT >::send(nest::Event& e, nest::thread t, double t_lastspike, const STDPSymHomCommonProperties &cp)
{
  // synapse STDP depressing/facilitation dynamics

	double t_spike = e.get_stamp().get_ms();

  
  // t_lastspike_ = 0 initially
  nest::Node* target = get_target( t );
  double dendritic_delay = nest::Time(nest::Time::step(get_delay())).get_ms();
    
  //get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque<mynest::histentry_sym>::iterator start;
  std::deque<mynest::histentry_sym>::iterator finish;
  //((mynest::Archiving_Node_Sym *)target_)->get_sym_history(t_lastspike - dendritic_delay, t_spike - dendritic_delay,&start, &finish);
  ((mynest::Archiving_Node_Sym *)target)->get_sym_history(t_lastspike, t_spike,&start, &finish);
  //weight change due to post-synaptic spikes since last pre-synaptic spike
  double_t dt;
  while (start != finish)
  {

	//dt = t_lastspike - (start->t_ + dendritic_delay);
	dt = t_lastspike - start->t_;
	//double old_time = start->t_ + dendritic_delay;
	double old_time = start->t_;
    ++start;


//    std::cout << Kexpt1_ << std::endl;

    // Update trace at the time of the postsynaptic spike
    // Calculate central = exp(-abs(dt/tau1))*cos(dt*pi/(tau1*2))^2
    double dt_tau1 = dt*cp.inv_tau_sym1_;

    double dt_pi_tau1 = M_PI*dt_tau1;
    double aux_expon_tau1 = std::exp(dt_tau1);
    double aux_cos_tau1 = std::cos(dt_pi_tau1);
    double aux_sin_tau1 = std::sin(dt_pi_tau1);

    double exp_t1 = Kexpt1_*aux_expon_tau1;

    // Calculate cos_2_t1 = Cos(dt*pi/tau1)
    double cos_2_t1 = aux_expon_tau1*(Kcos2t1_*aux_cos_tau1 - Ksin2t1_*aux_sin_tau1);

    double central = 0.5*(exp_t1 + cos_2_t1);

//    if (central<0.0){
//		std::cout << "Error: Central<0.0. dt=" << dt << " exp_t1=" << exp_t1 << " and cos_2_t1=" << cos_2_t1 << std::endl;
//		std::cout << "Original at time " << t_lastspike << " exp_t1=" << Kexpt1_ << " Kcos2t1_=" << Kcos2t1_ << " Ksin2t1_=" << Ksin2t1_ << std::endl;
//		std::cout << "New at time " << old_time << " aux_expon_tau1=" << aux_expon_tau1 << " aux_cos_tau1=" << aux_cos_tau1 << " aux_sin_tau1=" << aux_sin_tau1 << std::endl;
//	}

    // Calculate external = A*exp(-2*abs(dt)/tau2)*sin(dt*pi/(tau2*2))^2
    double dt_tau2 = dt*cp.inv_tau_sym2_;
//    std::cout << "dt: " << dt << " - inv_tau_sym2: " << cp.inv_tau_sym2_ << std::endl;
    double dt_pi_tau2 = M_PI*dt_tau2;
    double aux_expon_tau2 = std::exp(2*dt_tau2);
    double aux_cos_tau2 = std::cos(dt_pi_tau2);
    double aux_sin_tau2 = std::sin(dt_pi_tau2);

    double exp_t2 = Kexpt2_*aux_expon_tau2;

    // Calculate cos_2_t2 = Cos(dt*pi/tau2)
    double cos_2_t2 = aux_expon_tau2*(Kcos2t2_*aux_cos_tau2 + Ksin2t2_*aux_sin_tau2);

    double external = 0.5*(exp_t2 - cos_2_t2);

//    if (external<0.0){
//		std::cout << "Error: External<0.0. dt=" << dt << " exp_t2=" << exp_t2 << " and cos_2_t2=" << cos_2_t2 << std::endl;
//		std::cout << "Original at time " << t_lastspike << " exp_t2=" << Kexpt2_ << " Kcos2t2_=" << Kcos2t2_ << " Ksin2t2_=" << Ksin2t2_ << std::endl;
//		std::cout << "New at time " << old_time << " aux_expon_tau2=" << aux_expon_tau2 << " aux_cos_tau2=" << aux_cos_tau2 << " aux_sin_tau2=" << aux_sin_tau2 << std::endl;
//
//	}

//    std::cout << "Pre-post at time: " << t_lastspike + dt << " Central: " << central << " External: " << cp.sym_A_*external << " Alpha: " << cp.alpha_ << std::endl;

    //std::cout << "Pre-post traces. Central: " << central << " Expt1: " << exp_t1 << " Cos2t1: " << cos_2_t1 << " External: " << external << " Expt2: " << exp_t2 << " Cos2t2: " << cos_2_t2 << std::endl;

    double trace = central - cp.alpha_ * cp.sym_A_ * external;
    double old_weight = weight_;
    weight_ = apply_weight_change_(weight_, trace, cp);

    //std::cout << "Weight change (pre->post) from pre at " << t_lastspike << " to post at time " << t_lastspike-dt << ": " << weight_-old_weight << ". Current value: " << weight_ << std::endl;

//    std::cout << "Post-spike received at dt" << dt << ": Trace: " << trace << " - Kexpt1 - " << exp_t1 << " Kcos2t1 - " << cos_2_t1 <<
//    	    	" Kexpt2 - " << exp_t2 << " Kcos2t2 - " << cos_2_t2 << std::endl;

  }

  //depression due to the incoming pre-synaptic spike
  double central, external;
  //((mynest::Archiving_Node_Sym *)target_)->get_sym_K_value(t_spike - dendritic_delay, central, external);
  ((mynest::Archiving_Node_Sym *)target)->get_sym_K_value(t_spike, central, external);
//  std::cout << "Post-pre at time: " << t_spike - dendritic_delay << " Central: " << central << " External: " << cp.sym_A_*external << " Alpha: " << cp.alpha_ << std::endl;

  double old_weight = weight_;
  weight_ = apply_weight_change_(weight_, central - cp.alpha_ * cp.sym_A_ * external, cp);

  //std::cout << "Weight change (post->pre) from pre at " << t_spike-dendritic_delay << " to post accumulated: " << weight_-old_weight << ". Current value: " << weight_ << std::endl;
  //std::cout << "Weight change (post->pre) from pre at " << t_spike << " to post accumulated: " << weight_-old_weight << ". Current value: " << weight_ << std::endl;


  e.set_receiver(*target);
  e.set_weight(weight_);
  e.set_delay(get_delay_steps());
  e.set_rport(get_rport());
  e();

//  std::cout << "Updating synapsis from time " << t_lastspike << " to " << t_spike << ". Old values: " << "- Kexpt1 - " << Kexpt1_ << " Kcos2t1 - " << Kcos2t1_ << " Ksin2t1 - " << Ksin2t1_ << " Kexpt2 - " << Kexpt2_ << " Ksin2t2 - " << Ksin2t2_ << " Kcos2t2 - " << Kcos2t2_ << this << std::endl;

  dt = t_lastspike - t_spike;

  // Calculate central = exp(-abs(dt/tau1))*cos(dt*pi/(tau1*2))^2
  double dt_tau1 = dt*cp.inv_tau_sym1_;

//  std::cout << "dt_tau1-" << dt_tau1 << std::endl;

  double dt_pi_tau1 = M_PI*dt_tau1;
  double aux_expon_tau1 = std::exp(dt_tau1);
  double aux_cos_tau1 = std::cos(dt_pi_tau1);
  double aux_sin_tau1 = std::sin(dt_pi_tau1);

//  std::cout << " aux_expon_tau1-" << aux_expon_tau1 << " aux_cos_tau1-" << aux_cos_tau1 << " aux_sin_tau1-" << aux_sin_tau1 << std::endl;

  Kexpt1_ = Kexpt1_*aux_expon_tau1 + 1.0;

  // Calculate cos_2_t1 = Cos(dt*pi/tau1)
  Kcos2t1_ = aux_expon_tau1*(Kcos2t1_*aux_cos_tau1 - Ksin2t1_*aux_sin_tau1) + 1.0;
  Ksin2t1_ = aux_expon_tau1*(Ksin2t1_*aux_cos_tau1 + Kcos2t1_*aux_sin_tau1);

  // std::cout << "Kexpt1-" << Kexpt1_ << "Kcos2t1-" << Kcos2t1_ << "Ksin2t1-" << Ksin2t1_ << std::endl;

  // Calculate external = A*exp(-2*abs(dt)/tau2)*sin(dt*pi/(tau2*2))^2
  double dt_tau2 = dt*cp.inv_tau_sym2_;
  double dt_pi_tau2 = M_PI*dt_tau2;
  double aux_expon_tau2 = std::exp(2*dt_tau2);
  double aux_cos_tau2 = std::cos(dt_pi_tau2);
  double aux_sin_tau2 = std::sin(dt_pi_tau2);

  Kexpt2_ = Kexpt2_*aux_expon_tau2 + 1.0;

  // Calculate cos_2_t2 = Cos(dt*pi/tau2)
  Kcos2t2_ = aux_expon_tau2*(Kcos2t2_*aux_cos_tau2 + Ksin2t2_*aux_sin_tau2) + 1.0;
  Ksin2t2_ = aux_expon_tau2*(Ksin2t2_*aux_cos_tau2 + Kcos2t2_*aux_sin_tau2);

  //std::cout << "Updated presynaptic trace from " << t_lastspike << " to " << t_spike << ". Expt1: " << Kexpt1_ << " Cos2t1: " << Kcos2t1_ << "Sin2t1: " << Ksin2t1_ << " Expt2: " << Kexpt2_ << " Cos2t2: " << Kcos2t2_ << " Sin2t2: " << Ksin2t2_ << std::endl;

//  std::cout << "Pre-spike received at " << t_spike << ": Trace: " << central - cp.alpha_ * cp.sym_A_ * external << "- Kexpt1 - " << Kexpt1_ << " Kcos2t1 - " << Kcos2t1_ << " Ksin2t1 - " << Ksin2t1_ << " dt_tau1 - " << dt_tau1 << " dt_pi_tau1 - " << dt_pi_tau1 << " aux_expon_tau1 - " << aux_expon_tau1 << " aux_cos_tau1 - " << aux_cos_tau1 << " aux_sin_tau1 - " << aux_sin_tau1
//		  	  	  << " Kexpt2 - " << Kexpt2_ << " Ksin2t2 - " << Ksin2t2_ << " Kcos2t2 - " << Kcos2t2_ << " dt_tau2 - " << dt_tau2 << " dt_pi_tau2 - " << dt_pi_tau2 << " aux_expon_tau2 - " << aux_expon_tau2 << " aux_cos_tau2 - " << aux_cos_tau2 << " aux_sin_tau2 - " << aux_sin_tau2 << std::endl;
//  std::cout << "Calculating new STDPSym values: Kexpt1-" << Kexpt1_ << " Kcos2t1-" << Kcos2t1_ << " Ksin2t1-" << Ksin2t1_ << " Kexpt2-" << Kexpt2_ << " Ksin2t2-" << Ksin2t2_ << " Kcos2t2-" << Kcos2t2_ << std::endl;

  }

template < typename targetidentifierT >
void STDPSymConnectionHom< targetidentifierT >::get_status(DictionaryDatum & d) const
{

    // base class properties, different for individual synapse
    ConnectionBase::get_status(d);

    // own properties, different for individual synapse
    def<double>(d, "Kexpt1", Kexpt1_);
    def<double>(d, "Kcos2t1", Kcos2t1_);
    def<double>(d, "Ksin2t1", Ksin2t1_);
    def<double>(d, "Kexpt2", Kexpt2_);
    def<double>(d, "Kcos2t2", Kcos2t2_);
    def<double>(d, "Ksin2t2", Ksin2t2_);
    def<double>(d, nest::names::weight, weight_ );
    def<long>(d, nest::names::size_of, sizeof(*this));
  }

template < typename targetidentifierT >
void STDPSymConnectionHom< targetidentifierT >::set_status(const DictionaryDatum & d, nest::ConnectorModel &cm)
{
    // base class properties
    ConnectionBase::set_status(d, cm);
    updateValue<double>(d, "Kexpt1", Kexpt1_);
    updateValue<double>(d, "Kcos2t1", Kcos2t1_);
    updateValue<double>(d, "Ksin2t1", Ksin2t1_);
    updateValue<double>(d, "Kexpt2", Kexpt2_);
    updateValue<double>(d, "Kcos2t2", Kcos2t2_);
    updateValue<double>(d, "Ksin2t2", Ksin2t2_);
    updateValue<double>(d, nest::names::weight, weight_);

//    std::cout << "Updating STDPSym values: Kexpt1-" << Kexpt1_ << " Kcos2t1-" << Kcos2t1_ << " Ksin2t1-" << Ksin2t1_ << " Kexpt2-" << Kexpt2_ << " Ksin2t2-" << Ksin2t2_ << " Kcos2t2-" << Kcos2t2_ << this << std::endl;

    // exception throwing must happen after setting own parameters
    // this exception will be caught and ignored within generic_connector_model::set_status()
    // it will not be caught if set_status is called directly to signify the specified error
    //if (d->known("tau_sym") || d->known("lambda") || d->known("alpha") || d->known("Wmax") ){
	//  throw nest::ChangeCommonPropsByIndividual("STDPSymConnectionHom::set_status(): you are trying to set common properties via an individual synapse.");
    //}

  }

} // of namespace nest

#endif // of #ifndef STDP_CONNECTION_HOM_H
