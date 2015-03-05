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

#include "connection_het_wd.h"
#include "archiving_node_sym.h"
#include <cmath>

namespace mynest
{

  /**
   * Class containing the common properties for all synapses of type STDPConnectionHom.
   */
  class STDPSymHomCommonProperties : public nest::CommonSynapseProperties
    {
      friend class STDPSymConnectionHom;

    public:

      /**
       * Default constructor.
       * Sets all property values to defaults.
       */
      STDPSymHomCommonProperties();
   
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
      nest::double_t inv_tau_sym1_;
      nest::double_t lambda_;
      nest::double_t alpha_;
      nest::double_t Wmax_;

      // Auxiliar variables common to all connections
      nest::double_t inv_tau_sym2_;
      nest::double_t sym_A_;
    };



  /**
   * Class representing an STDP connection with homogeneous parameters, i.e. parameters are the same for all synapses.
   */
  class STDPSymConnectionHom : public nest::ConnectionHetWD
  {

  public:
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

  /**
   * Default Destructor.
   */
  virtual ~STDPSymConnectionHom() {}

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
  void send(nest::Event& e, nest::double_t t_lastspike, const STDPSymHomCommonProperties &);

  // overloaded for all supported event types
  using Connection::check_event;
  void check_event(nest::SpikeEvent&) {}
  
 private:

  nest::double_t apply_weight_change_(nest::double_t w, nest::double_t trace, const STDPSymHomCommonProperties &cp);

  // data members of each connection
  // Accumulation variables
  nest::double_t Kexpt1_; 		// value of exp(-abs(x)/tau1) at that time
  nest::double_t Kcos2t1_;		// value of exp(-abs(x)/tau1)*cos(2*x/tau1) at that time
  nest::double_t Ksin2t1_;		// value of exp(-abs(x)/tau1)*sin(2*x/tau1) at that time
  nest::double_t Kexpt2_; 		// value of exp(-abs(x)/tau1) at that time
  nest::double_t Ksin2t2_;		// value of exp(-abs(x)/tau2)*sin(2*x/tau2) at that time
  nest::double_t Kcos2t2_;		// value of exp(-abs(x)/tau2)*cos(2*x/tau2) at that time

  };


inline
nest::double_t STDPSymConnectionHom::apply_weight_change_(nest::double_t w, nest::double_t trace, const STDPSymHomCommonProperties &cp)
{
  nest::double_t norm_w = (w / cp.Wmax_) + (cp.lambda_ * trace);
//  std::cout << "Old weight: " << w << " Norm. Weight: " << norm_w << " Kexpt1: " << Kexpt1_ << " Kcos2t1: " << Kcos2t1_ << " Ksin2t1: " << Ksin2t1_ << " Kexpt2: " << Kexpt2_ << " Ksin2t2: " << Ksin2t2_ << " Kcos2t2: " << Kcos2t2_ << " Trace: " << trace << " Lambda: " << cp.lambda_ << " Wmax: " << cp.Wmax_ << this << std::endl;
    if (norm_w > 1.0){
	  return cp.Wmax_;
    }else if (norm_w < 0.0)
	  return 0.0;

  return norm_w * cp.Wmax_;
}

inline 
  void STDPSymConnectionHom::check_connection(nest::Node & s, nest::Node & r, nest::rport receptor_type, nest::double_t t_lastspike)
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
void STDPSymConnectionHom::send(nest::Event& e, nest::double_t t_lastspike, const STDPSymHomCommonProperties &cp)
{
  // synapse STDP depressing/facilitation dynamics

	nest::double_t t_spike = e.get_stamp().get_ms();

  
  // t_lastspike_ = 0 initially
  
  //std::cout << "Spike sent to " << target_ << " in time " << t_spike << std::endl;

	nest::double_t dendritic_delay = nest::Time(nest::Time::step(delay_)).get_ms();
    
  //get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque<mynest::histentry_sym>::iterator start;
  std::deque<mynest::histentry_sym>::iterator finish;
  //((mynest::Archiving_Node_Sym *)target_)->get_sym_history(t_lastspike - dendritic_delay, t_spike - dendritic_delay,&start, &finish);
  ((mynest::Archiving_Node_Sym *)target_)->get_sym_history(t_lastspike, t_spike,&start, &finish);
  //weight change due to post-synaptic spikes since last pre-synaptic spike
  double_t dt;
  while (start != finish)
  {

	//dt = t_lastspike - (start->t_ + dendritic_delay);
	dt = t_lastspike - start->t_;
	//nest::double_t old_time = start->t_ + dendritic_delay;
	nest::double_t old_time = start->t_;
    ++start;


//    std::cout << Kexpt1_ << std::endl;

    // Update trace at the time of the postsynaptic spike
    // Calculate central = exp(-abs(dt/tau1))*cos(dt*pi/(tau1*2))^2
    nest::double_t dt_tau1 = dt*cp.inv_tau_sym1_;

    nest::double_t dt_pi_tau1 = M_PI*dt_tau1;
    nest::double_t aux_expon_tau1 = std::exp(dt_tau1);
    nest::double_t aux_cos_tau1 = std::cos(dt_pi_tau1);
    nest::double_t aux_sin_tau1 = std::sin(dt_pi_tau1);

    nest::double_t exp_t1 = Kexpt1_*aux_expon_tau1;

    // Calculate cos_2_t1 = Cos(dt*pi/tau1)
    nest::double_t cos_2_t1 = aux_expon_tau1*(Kcos2t1_*aux_cos_tau1 - Ksin2t1_*aux_sin_tau1);

    nest::double_t central = 0.5*(exp_t1 + cos_2_t1);

//    if (central<0.0){
//		std::cout << "Error: Central<0.0. dt=" << dt << " exp_t1=" << exp_t1 << " and cos_2_t1=" << cos_2_t1 << std::endl;
//		std::cout << "Original at time " << t_lastspike << " exp_t1=" << Kexpt1_ << " Kcos2t1_=" << Kcos2t1_ << " Ksin2t1_=" << Ksin2t1_ << std::endl;
//		std::cout << "New at time " << old_time << " aux_expon_tau1=" << aux_expon_tau1 << " aux_cos_tau1=" << aux_cos_tau1 << " aux_sin_tau1=" << aux_sin_tau1 << std::endl;
//	}

    // Calculate external = A*exp(-2*abs(dt)/tau2)*sin(dt*pi/(tau2*2))^2
    nest::double_t dt_tau2 = dt*cp.inv_tau_sym2_;
//    std::cout << "dt: " << dt << " - inv_tau_sym2: " << cp.inv_tau_sym2_ << std::endl;
    nest::double_t dt_pi_tau2 = M_PI*dt_tau2;
    nest::double_t aux_expon_tau2 = std::exp(2*dt_tau2);
    nest::double_t aux_cos_tau2 = std::cos(dt_pi_tau2);
    nest::double_t aux_sin_tau2 = std::sin(dt_pi_tau2);

    nest::double_t exp_t2 = Kexpt2_*aux_expon_tau2;

    // Calculate cos_2_t2 = Cos(dt*pi/tau2)
    nest::double_t cos_2_t2 = aux_expon_tau2*(Kcos2t2_*aux_cos_tau2 + Ksin2t2_*aux_sin_tau2);

    nest::double_t external = 0.5*(exp_t2 - cos_2_t2);

//    if (external<0.0){
//		std::cout << "Error: External<0.0. dt=" << dt << " exp_t2=" << exp_t2 << " and cos_2_t2=" << cos_2_t2 << std::endl;
//		std::cout << "Original at time " << t_lastspike << " exp_t2=" << Kexpt2_ << " Kcos2t2_=" << Kcos2t2_ << " Ksin2t2_=" << Ksin2t2_ << std::endl;
//		std::cout << "New at time " << old_time << " aux_expon_tau2=" << aux_expon_tau2 << " aux_cos_tau2=" << aux_cos_tau2 << " aux_sin_tau2=" << aux_sin_tau2 << std::endl;
//
//	}

//    std::cout << "Pre-post at time: " << t_lastspike + dt << " Central: " << central << " External: " << cp.sym_A_*external << " Alpha: " << cp.alpha_ << std::endl;

    //std::cout << "Pre-post traces. Central: " << central << " Expt1: " << exp_t1 << " Cos2t1: " << cos_2_t1 << " External: " << external << " Expt2: " << exp_t2 << " Cos2t2: " << cos_2_t2 << std::endl;

    nest::double_t trace = central - cp.alpha_ * cp.sym_A_ * external;
    nest::double_t old_weight = weight_;
    weight_ = apply_weight_change_(weight_, trace, cp);

    //std::cout << "Weight change (pre->post) from pre at " << t_lastspike << " to post at time " << t_lastspike-dt << ": " << weight_-old_weight << ". Current value: " << weight_ << std::endl;

//    std::cout << "Post-spike received at dt" << dt << ": Trace: " << trace << " - Kexpt1 - " << exp_t1 << " Kcos2t1 - " << cos_2_t1 <<
//    	    	" Kexpt2 - " << exp_t2 << " Kcos2t2 - " << cos_2_t2 << std::endl;

  }

  //depression due to the incoming pre-synaptic spike
  nest::double_t central, external;
  //((mynest::Archiving_Node_Sym *)target_)->get_sym_K_value(t_spike - dendritic_delay, central, external);
  ((mynest::Archiving_Node_Sym *)target_)->get_sym_K_value(t_spike, central, external);
//  std::cout << "Post-pre at time: " << t_spike - dendritic_delay << " Central: " << central << " External: " << cp.sym_A_*external << " Alpha: " << cp.alpha_ << std::endl;

  nest::double_t old_weight = weight_;
  weight_ = apply_weight_change_(weight_, central - cp.alpha_ * cp.sym_A_ * external, cp);

  //std::cout << "Weight change (post->pre) from pre at " << t_spike-dendritic_delay << " to post accumulated: " << weight_-old_weight << ". Current value: " << weight_ << std::endl;
  //std::cout << "Weight change (post->pre) from pre at " << t_spike << " to post accumulated: " << weight_-old_weight << ". Current value: " << weight_ << std::endl;


  e.set_receiver(*target_);
  e.set_weight(weight_);
  e.set_delay(delay_);
  e.set_rport(rport_);
  e();

//  std::cout << "Updating synapsis from time " << t_lastspike << " to " << t_spike << ". Old values: " << "- Kexpt1 - " << Kexpt1_ << " Kcos2t1 - " << Kcos2t1_ << " Ksin2t1 - " << Ksin2t1_ << " Kexpt2 - " << Kexpt2_ << " Ksin2t2 - " << Ksin2t2_ << " Kcos2t2 - " << Kcos2t2_ << this << std::endl;

  dt = t_lastspike - t_spike;

  // Calculate central = exp(-abs(dt/tau1))*cos(dt*pi/(tau1*2))^2
  nest::double_t dt_tau1 = dt*cp.inv_tau_sym1_;

//  std::cout << "dt_tau1-" << dt_tau1 << std::endl;

  nest::double_t dt_pi_tau1 = M_PI*dt_tau1;
  nest::double_t aux_expon_tau1 = std::exp(dt_tau1);
  nest::double_t aux_cos_tau1 = std::cos(dt_pi_tau1);
  nest::double_t aux_sin_tau1 = std::sin(dt_pi_tau1);

//  std::cout << " aux_expon_tau1-" << aux_expon_tau1 << " aux_cos_tau1-" << aux_cos_tau1 << " aux_sin_tau1-" << aux_sin_tau1 << std::endl;

  Kexpt1_ = Kexpt1_*aux_expon_tau1 + 1.0;

  // Calculate cos_2_t1 = Cos(dt*pi/tau1)
  Kcos2t1_ = aux_expon_tau1*(Kcos2t1_*aux_cos_tau1 - Ksin2t1_*aux_sin_tau1) + 1.0;
  Ksin2t1_ = aux_expon_tau1*(Ksin2t1_*aux_cos_tau1 + Kcos2t1_*aux_sin_tau1);

  // std::cout << "Kexpt1-" << Kexpt1_ << "Kcos2t1-" << Kcos2t1_ << "Ksin2t1-" << Ksin2t1_ << std::endl;

  // Calculate external = A*exp(-2*abs(dt)/tau2)*sin(dt*pi/(tau2*2))^2
  nest::double_t dt_tau2 = dt*cp.inv_tau_sym2_;
  nest::double_t dt_pi_tau2 = M_PI*dt_tau2;
  nest::double_t aux_expon_tau2 = std::exp(2*dt_tau2);
  nest::double_t aux_cos_tau2 = std::cos(dt_pi_tau2);
  nest::double_t aux_sin_tau2 = std::sin(dt_pi_tau2);

  Kexpt2_ = Kexpt2_*aux_expon_tau2 + 1.0;

  // Calculate cos_2_t2 = Cos(dt*pi/tau2)
  Kcos2t2_ = aux_expon_tau2*(Kcos2t2_*aux_cos_tau2 + Ksin2t2_*aux_sin_tau2) + 1.0;
  Ksin2t2_ = aux_expon_tau2*(Ksin2t2_*aux_cos_tau2 + Kcos2t2_*aux_sin_tau2);

  //std::cout << "Updated presynaptic trace from " << t_lastspike << " to " << t_spike << ". Expt1: " << Kexpt1_ << " Cos2t1: " << Kcos2t1_ << "Sin2t1: " << Ksin2t1_ << " Expt2: " << Kexpt2_ << " Cos2t2: " << Kcos2t2_ << " Sin2t2: " << Ksin2t2_ << std::endl;

//  std::cout << "Pre-spike received at " << t_spike << ": Trace: " << central - cp.alpha_ * cp.sym_A_ * external << "- Kexpt1 - " << Kexpt1_ << " Kcos2t1 - " << Kcos2t1_ << " Ksin2t1 - " << Ksin2t1_ << " dt_tau1 - " << dt_tau1 << " dt_pi_tau1 - " << dt_pi_tau1 << " aux_expon_tau1 - " << aux_expon_tau1 << " aux_cos_tau1 - " << aux_cos_tau1 << " aux_sin_tau1 - " << aux_sin_tau1
//		  	  	  << " Kexpt2 - " << Kexpt2_ << " Ksin2t2 - " << Ksin2t2_ << " Kcos2t2 - " << Kcos2t2_ << " dt_tau2 - " << dt_tau2 << " dt_pi_tau2 - " << dt_pi_tau2 << " aux_expon_tau2 - " << aux_expon_tau2 << " aux_cos_tau2 - " << aux_cos_tau2 << " aux_sin_tau2 - " << aux_sin_tau2 << std::endl;
//  std::cout << "Calculating new STDPSym values: Kexpt1-" << Kexpt1_ << " Kcos2t1-" << Kcos2t1_ << " Ksin2t1-" << Ksin2t1_ << " Kexpt2-" << Kexpt2_ << " Ksin2t2-" << Ksin2t2_ << " Kcos2t2-" << Kcos2t2_ << std::endl;

  }

} // of namespace nest

#endif // of #ifndef STDP_CONNECTION_HOM_H
