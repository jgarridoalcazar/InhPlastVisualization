/*
 *  iaf_cond_exp_ip_sym.h
 *
 *  This file is based on the iaf_cond_exp cell model distributed with NEST.
 *  
 *  Modified by: Jes�s Garrido (jgarridoalcazar at gmail.com) in 2014.
 */

#ifndef IAF_COND_EXP_IP_SYM_H
#define IAF_COND_EXP_IP_SYM_H

#include "config.h"

#ifdef HAVE_GSL

#include "nest.h"
#include "event.h"
#include "archiving_node_sym.h"
#include "ring_buffer.h"
#include "connection.h"
#include "universal_data_logger.h"
#include "recordables_map.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

/* BeginDocumentation
Name: iaf_cond_exp_ip_sym - Conductance based leaky integrate-and-fire neuron model with intrinsic plasticity.

Description:
iaf_cond_exp_sym is an implementation of a spiking neuron using IAF dynamics with
conductance-based synapses and intrinsic plasticity. Incoming spike events induce a post-synaptic change 
of conductance modelled by an exponential function. The exponential function 
is normalised such that an event of weight 1.0 results in a peak conductance of 1 nS.
The intrinsic plasticity mechanism adjust the membrane capacitance and resistance in order to keep
constant firing rate and maximizing the sensitivity to the input conductances. This plasticity
is based on the paper [1]. 

Parameters: 
The following parameters can be set in the status dictionary.

V_m        double - Membrane potential in mV 
E_L        double - Leak reversal potential in mV.
r_C         double - Inverse of the capacity of the membrane in pF-1
t_ref      double - Duration of refractory period in ms. 
V_th       double - Spike threshold in mV.
V_reset    double - Reset potential of the membrane in mV.
E_ex       double - Excitatory reversal potential in mV.
E_in       double - Inhibitory reversal potential in mV.
g_L        double - Leak conductance in nS (also known rR in [1]);
tau_syn_ex double - Time constant of the excitatory synaptic exponential function in ms.
tau_syn_in double - Time constant of the inhibitory synaptic exponential function in ms.
I_e        double - Constant external input current in pA.
tau_ip	   double - Intrinsic plasticity time constant in ms
epsilon_rC double - Effect of each postsynaptic spike in the rC state variable
epsilon_rR double - Effect of each postsynaptic spike in the rR state variable
beta       double - Parameter of the probability Weibull distribution.

Sends: SpikeEvent

Receives: SpikeEvent, CurrentEvent, DataLoggingRequest

References: 

[1] Li, C., & Li, Y. (2013). A spike-based model of neuronal intrinsic plasticity.
Autonomous Mental Development, IEEE Transactions on, 5(1), 62-73.

Author: Jesus Garrido

SeeAlso: iaf_cond_exp
*/

// Define name constants for state variables and parameters
namespace nest
{
	namespace names
	{
    	// Neuron parameters
    	extern const Name r_C;                 	//!< Inverse of the capacity of the membrane
    	extern const Name tau_ip;				//!< IP time constant
    	extern const Name epsilon_rC;			//!< IP effect of postsynaptic spike in the rC variable
    	extern const Name epsilon_rR;			//!< IP effect of postsynaptic spike in the g_L variable
    	extern const Name beta;					//!< IP parameter of the probability distribution
	}
}

namespace mynest
{
  /**
   * Function computing right-hand side of ODE for GSL solver.
   * @note Must be declared here so we can befriend it in class.
   * @note Must have C-linkage for passing to GSL. Internally, it is
   *       a first-class C++ function, but cannot be a member function
   *       because of the C-linkage.
   * @note No point in declaring it inline, since it is called
   *       through a function pointer.
   * @param void* Pointer to model neuron instance.
   */
  extern "C"
  int iaf_cond_exp_ip_sym_dynamics (double, const double*, double*, void*);
  
  class iaf_cond_exp_ip_sym : public mynest::Archiving_Node_Sym
  {
    
  public:        
    
    iaf_cond_exp_ip_sym();
    iaf_cond_exp_ip_sym(const iaf_cond_exp_ip_sym&);
    ~iaf_cond_exp_ip_sym();

    /**
     * Import sets of overloaded virtual functions.
     * We need to explicitly include sets of overloaded
     * virtual functions into the current scope.
     * According to the SUN C++ FAQ, this is the correct
     * way of doing things, although all other compilers
     * happily live without.
     */

    using nest::Node::handles_test_event;
    using nest::Node::handle;

    nest::port send_test_event(nest::Node&, nest::rport, nest::synindex, bool);
    
    void handle(nest::SpikeEvent &);
    void handle(nest::CurrentEvent &);
    void handle(nest::DataLoggingRequest &); 
    
    nest::port handles_test_event(nest::SpikeEvent &, nest::rport);
    nest::port handles_test_event(nest::CurrentEvent &, nest::rport);
    nest::port handles_test_event(nest::DataLoggingRequest &, nest::rport);
    
    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);
    
  private:
    void init_state_(const Node& proto);
    void init_buffers_();
    void calibrate();
    void update(nest::Time const &, const nest::long_t, const nest::long_t);

    // END Boilerplate function declarations ----------------------------

    // Friends --------------------------------------------------------

    // make dynamics function quasi-member
    friend int iaf_cond_exp_ip_sym_dynamics(double, const double*, double*, void*);

    // The next two classes need to be friends to access the State_ class/member
    friend class nest::RecordablesMap<iaf_cond_exp_ip_sym>;
    friend class nest::UniversalDataLogger<iaf_cond_exp_ip_sym>;

  private:

    // ---------------------------------------------------------------- 

    //! Model parameters
    struct Parameters_ {
      nest::double_t V_th_;       //!< Threshold Potential in mV
      nest::double_t V_reset_;    //!< Reset Potential in mV
      nest::double_t t_ref_;      //!< Refractory period in ms
      nest::double_t g_L;         //!< Initial leak Conductance in nS (also known as rR)
      nest::double_t r_C;         //!< Initial inverse of the membrane Capacitance in pF
      nest::double_t E_ex;        //!< Excitatory reversal Potential in mV
      nest::double_t E_in;        //!< Inhibitory reversal Potential in mV
      nest::double_t E_L;         //!< Leak reversal Potential (aka resting potential) in mV
      nest::double_t tau_synE;    //!< Synaptic Time Constant Excitatory Synapse in ms
      nest::double_t tau_synI;    //!< Synaptic Time Constant for Inhibitory Synapse in ms
      nest::double_t I_e;         //!< Constant Current in pA
      nest::double_t tau_ip;      //!< Intrinsic plasticity time constant in ms
      nest::double_t epsilon_rC;	//!< Effect of each postsynaptic spike in the rC state variable
      nest::double_t epsilon_rR;	//!< Effect of each postsynaptic spike in the g_L state variable
      nest::double_t beta;		//!< Parameter of the probability Weibull distribution describing the input current.
    
      Parameters_();  //!< Sets default parameter values

      void get(DictionaryDatum&) const;  //!< Store current values in dictionary
      void set(const DictionaryDatum&);  //!< Set values from dicitonary
    };

  public:
    // ---------------------------------------------------------------- 

    /**
     * State variables of the model.
     * @note Copy constructor and assignment operator required because
     *       of C-style array.
     */
    struct State_ {

      //! Symbolic indices to the elements of the state vector y
      enum StateVecElems { V_M = 0,           
			   G_EXC,     
			   G_INH,
			   G_L,
			   R_C,
			   STATE_VEC_SIZE };

      nest::double_t y_[STATE_VEC_SIZE];  //!< neuron state, must be C-array for GSL solver
      nest::int_t    r_;                  //!< number of refractory steps remaining

      State_(const Parameters_&);  //!< Default initialization
      State_(const State_&);
      State_& operator=(const State_&);

      void get(DictionaryDatum&) const;
      void set(const DictionaryDatum&, const Parameters_&);
    };    

    // ---------------------------------------------------------------- 

  private:
    /**
     * Buffers of the model.
     */
    struct Buffers_ {
      Buffers_(iaf_cond_exp_ip_sym&);                   //!<Sets buffer pointers to 0
      Buffers_(const Buffers_&, iaf_cond_exp_ip_sym&);  //!<Sets buffer pointers to 0

      //! Logger for all analog data
      nest::UniversalDataLogger<iaf_cond_exp_ip_sym> logger_;

      /** buffers and sums up incoming spikes/currents */
      nest::RingBuffer spike_exc_;
      nest::RingBuffer spike_inh_;
      nest::RingBuffer currents_;

      /** GSL ODE stuff */
      gsl_odeiv_step*    s_;    //!< stepping function
      gsl_odeiv_control* c_;    //!< adaptive stepsize control function
      gsl_odeiv_evolve*  e_;    //!< evolution function
      gsl_odeiv_system   sys_;  //!< struct describing system
      
      // IntergrationStep_ should be reset with the neuron on ResetNetwork,
      // but remain unchanged during calibration. Since it is initialized with
      // step_, and the resolution cannot change after nodes have been created,
      // it is safe to place both here.
      nest::double_t step_;           //!< step size in ms
      double   IntegrationStep_;//!< current integration time step, updated by GSL

      /** 
       * Input current injected by CurrentEvent.
       * This variable is used to transport the current applied into the
       * _dynamics function computing the derivative of the state vector.
       * It must be a part of Buffers_, since it is initialized once before
       * the first simulation, but not modified before later Simulate calls.
       */
      nest::double_t I_stim_;
    };

     // ---------------------------------------------------------------- 

     /**
      * Internal variables of the model.
      */
     struct Variables_ { 
    	nest::int_t    RefractoryCounts_;
     };

    // Access functions for UniversalDataLogger -------------------------------
    
    //! Read out state vector elements, used by UniversalDataLogger
    template <State_::StateVecElems elem>
    nest::double_t get_y_elem_() const { return S_.y_[elem]; }

    // ---------------------------------------------------------------- 

    Parameters_ P_;
    State_      S_;
    Variables_  V_;
    Buffers_    B_;

    //! Mapping of recordables names to access functions
    static nest::RecordablesMap<iaf_cond_exp_ip_sym> recordablesMap_;
  };

  
  inline nest::port iaf_cond_exp_ip_sym::send_test_event(nest::Node& target, nest::rport receptor_type, nest::synindex, bool)
  {
    nest::SpikeEvent e;
    e.set_sender(*this);
    return target.handles_test_event(e, receptor_type);
  }

  inline
  nest::port iaf_cond_exp_ip_sym::handles_test_event(nest::SpikeEvent&, nest::rport receptor_type)
  {
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return 0;
  }
 
  inline
  nest::port iaf_cond_exp_ip_sym::handles_test_event(nest::CurrentEvent&, nest::rport receptor_type)
  {
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return 0;
  }

  inline
  nest::port iaf_cond_exp_ip_sym::handles_test_event(nest::DataLoggingRequest& dlr,
		  nest::rport receptor_type)
  {
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return B_.logger_.connect_logging_device(dlr, recordablesMap_);
  }
 
  inline
  void iaf_cond_exp_ip_sym::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    S_.get(d);
    mynest::Archiving_Node_Sym::get_status(d);

    (*d)[nest::names::recordables] = recordablesMap_.get_list();
  }

  inline
  void iaf_cond_exp_ip_sym::set_status(const DictionaryDatum &d)
  {
    Parameters_ ptmp = P_;  // temporary copy in case of errors
    ptmp.set(d);                       // throws if BadProperty
    State_      stmp = S_;  // temporary copy in case of errors
    stmp.set(d, ptmp);                 // throws if BadProperty

    // We now know that (ptmp, stmp) are consistent. We do not 
    // write them back to (P_, S_) before we are also sure that 
    // the properties to be set in the parent class are internally 
    // consistent.
    mynest::Archiving_Node_Sym::set_status(d);

    // if we get here, temporaries contain consistent set of properties
    P_ = ptmp;
    S_ = stmp;
  }
  
} // namespace

#endif //HAVE_GSL
#endif //IAF_COND_EXP_H
