/*
 *  iaf_cond_exp_sto_ip.cpp
 *
 *  This file is based on the iaf_cond_exp cell model distributed with NEST.
 *  
 *  Modified by: Jes�s Garrido (jgarridoalcazar at gmail.com) in 2014.
 */

#include "iaf_cond_exp_sto_ip.h"

#ifdef HAVE_GSL

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include <limits>

#include "universal_data_logger_impl.h"
#include "event.h"

#include <iomanip>
#include <iostream>
#include <cstdio>

/* ---------------------------------------------------------------- 
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap<mynest::iaf_cond_exp_sto_ip> mynest::iaf_cond_exp_sto_ip::recordablesMap_;


namespace nest  // template specialization must be placed in namespace
{
  // Override the create() method with one call to RecordablesMap::insert_() 
  // for each quantity to be recorded.
  template <>
  void RecordablesMap<mynest::iaf_cond_exp_sto_ip>::create()
  {
    // use standard names whereever you can for consistency!
    insert_(names::V_m, 
	    &mynest::iaf_cond_exp_sto_ip::get_y_elem_<mynest::iaf_cond_exp_sto_ip::State_::V_M>);
    insert_(names::g_ex, 
	    &mynest::iaf_cond_exp_sto_ip::get_y_elem_<mynest::iaf_cond_exp_sto_ip::State_::G_EXC>);
    insert_(names::g_in, 
	    &mynest::iaf_cond_exp_sto_ip::get_y_elem_<mynest::iaf_cond_exp_sto_ip::State_::G_INH>);
    insert_(names::V_th,
    	    &mynest::iaf_cond_exp_sto_ip::get_y_elem_<mynest::iaf_cond_exp_sto_ip::State_::V_TH>);
    insert_(names::r_0,
    	    &mynest::iaf_cond_exp_sto_ip::get_y_elem_<mynest::iaf_cond_exp_sto_ip::State_::R_0>);
    insert_(names::u_alpha,
    	    &mynest::iaf_cond_exp_sto_ip::get_y_elem_<mynest::iaf_cond_exp_sto_ip::State_::U_ALPHA>);
    insert_(names::refractoriness,
        	    &mynest::iaf_cond_exp_sto_ip::get_y_elem_<mynest::iaf_cond_exp_sto_ip::State_::REFR>);
    insert_(names::gain,
        	    &mynest::iaf_cond_exp_sto_ip::get_y_elem_<mynest::iaf_cond_exp_sto_ip::State_::GAIN>);
    insert_(names::firing_probability,
        	    &mynest::iaf_cond_exp_sto_ip::get_y_elem_<mynest::iaf_cond_exp_sto_ip::State_::FIR_PROB>);
  }
  
  namespace names
  {

  	  const Name r_0("r_0");
      const Name u_alpha("u_alpha");
      const Name ip_rate("ip_rate");
      const Name target_firing("target_firing");

      // Debugging names
      const Name refractoriness("refractoriness");
      const Name gain("gain");
      const Name firing_probability("firing_probability");
  }
}

extern "C"
inline int mynest::iaf_cond_exp_sto_ip_dynamics(double, const double y[], double f[], void* pnode)
{ 
  // a shorthand
  typedef mynest::iaf_cond_exp_sto_ip::State_ S;

  // get access to node so we can almost work as in a member function
  assert(pnode);
  const mynest::iaf_cond_exp_sto_ip& node =  *(reinterpret_cast<mynest::iaf_cond_exp_sto_ip*>(pnode));

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[]. 

  // The following code is verbose for the sake of clarity. We assume that a
  // good compiler will optimize the verbosity away ...
  const double I_syn_exc = y[S::G_EXC] * (y[S::V_M] - node.P_.E_ex); 
  const double I_syn_inh = y[S::G_INH] * (y[S::V_M] - node.P_.E_in); 
  const double I_L       = node.P_.g_L * ( y[S::V_M] - node.P_.E_L );
  const double I_total   = node.P_.I_e - I_syn_exc - I_syn_inh;
  
  //V dot
  f[S::V_M] = ( - I_L + node.B_.I_stim_ + I_total) / node.P_.C_m; // Vm diff. equation
  f[S::G_EXC] = -y[S::G_EXC] / node.P_.tau_synE; // Gexc diff. equation
  f[S::G_INH] = -y[S::G_INH] / node.P_.tau_synI; // Ginh diff. equation
  f[S::V_TH] = f[S::R_0] = f[S::U_ALPHA] = f[S::REFR] = f[S::GAIN] = f[S::FIR_PROB] = 0.0;

  return GSL_SUCCESS;
}

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
mynest::iaf_cond_exp_sto_ip::Parameters_::Parameters_()
  : V_reset_   (-70.0    ),  // mV
    t_ref_abs  (  3.0    ),  // ms
    t_ref	   ( 10.0    ),  // ms
    g_L        ( 16.6667 ),  // nS
    C_m        (  2.5    ),  // pF-1
    E_ex       (  0.0    ),  // mV
    E_in       (-85.0    ),  // mV
    E_L        (-70.0    ),  // mV
    tau_synE   (  0.2    ),  // ms
    tau_synI   (  2.0    ),  // ms
    I_e        (  0.0    ),  // pA
	ip_rate    (  1.0e-4 ),  // Unitless
	target_firing(10.0   )	 // Hz
{
}

mynest::iaf_cond_exp_sto_ip::State_::State_(const Parameters_& p)
  : r_(0), time_rel_ref(0)
{
  y_[V_M] = p.V_reset_;
  y_[G_EXC] = y_[G_INH] = 0;
  y_[V_TH] = -50.0;
  y_[R_0] = 2.0;
  y_[U_ALPHA] = 2.0;
  y_[REFR] = 0.;
  y_[GAIN] = 0.;
  y_[FIR_PROB] = 0.;
}

mynest::iaf_cond_exp_sto_ip::State_::State_(const State_& s)
  : r_(s.r_), time_rel_ref(s.time_rel_ref)
{
  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y_[i] = s.y_[i];
}

mynest::iaf_cond_exp_sto_ip::State_& mynest::iaf_cond_exp_sto_ip::State_::operator=(const State_& s)
{
  assert(this != &s);  // would be bad logical error in program
  
  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y_[i] = s.y_[i];
  r_ = s.r_;
  time_rel_ref = s.time_rel_ref;
  return *this;
}

/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void mynest::iaf_cond_exp_sto_ip::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d,nest::names::V_reset,      V_reset_);
  def<double>(d,nest::names::t_ref_abs,    t_ref_abs);
  def<double>(d,nest::names::t_ref,        t_ref);
  def<double>(d,nest::names::g_L,          g_L);
  def<double>(d,nest::names::E_L,          E_L); 
  def<double>(d,nest::names::E_ex,         E_ex);
  def<double>(d,nest::names::E_in,         E_in);
  def<double>(d,nest::names::C_m,          C_m);
  def<double>(d,nest::names::tau_syn_ex,   tau_synE);
  def<double>(d,nest::names::tau_syn_in,   tau_synI);
  def<double>(d,nest::names::I_e,          I_e);
  def<double>(d,nest::names::ip_rate,      ip_rate);
  def<double>(d,nest::names::target_firing,target_firing);
}

void mynest::iaf_cond_exp_sto_ip::Parameters_::set(const DictionaryDatum& d)
{
  // allow setting the membrane potential
  updateValue<double>(d,nest::names::V_reset, 	  V_reset_);
  updateValue<double>(d,nest::names::t_ref_abs,   t_ref_abs);
  updateValue<double>(d,nest::names::t_ref,       t_ref);
  updateValue<double>(d,nest::names::E_L,         E_L);
  
  updateValue<double>(d,nest::names::E_ex,        E_ex);
  updateValue<double>(d,nest::names::E_in,        E_in);

  updateValue<double>(d,nest::names::C_m,         C_m);
  updateValue<double>(d,nest::names::g_L,         g_L);

  updateValue<double>(d,nest::names::tau_syn_ex, tau_synE);
  updateValue<double>(d,nest::names::tau_syn_in, tau_synI);

  updateValue<double>(d,nest::names::I_e,        I_e);

  updateValue<double>(d,nest::names::ip_rate,    ip_rate);
  updateValue<double>(d,nest::names::target_firing,target_firing);
  
  if ( C_m <= 0 )
    throw nest::BadProperty("Membrane capacitance must be strictly positive.");
    
  if ( t_ref < 0 )
    throw nest::BadProperty("Relative refractory time cannot be negative.");

  if ( t_ref_abs < 0 )
      throw nest::BadProperty("Absolute refractory time cannot be negative.");

  if ( tau_synE <= 0 || tau_synI <= 0)
    throw nest::BadProperty("All synpatic time constants must be strictly positive.");

  if (ip_rate < 0)
	throw nest::BadProperty("IP time constant must be positive or zero.");

  if (target_firing <= 0)
	throw nest::BadProperty("Target firing frequency must be positive.");
}

void mynest::iaf_cond_exp_sto_ip::State_::get(DictionaryDatum &d) const
{
  def<double>(d, nest::names::V_m, y_[V_M]); // Membrane potential
  def<double>(d, nest::names::V_th, y_[V_TH]);
  def<double>(d, nest::names::r_0, y_[R_0]);
  def<double>(d, nest::names::u_alpha, y_[U_ALPHA]);
}

void mynest::iaf_cond_exp_sto_ip::State_::set(const DictionaryDatum& d, const Parameters_&)
{
  updateValue<double>(d, nest::names::V_m, y_[V_M]);
  updateValue<double>(d, nest::names::V_th, y_[V_TH]);
  updateValue<double>(d, nest::names::r_0, y_[R_0]);
  updateValue<double>(d, nest::names::u_alpha, y_[U_ALPHA]);
}

mynest::iaf_cond_exp_sto_ip::Buffers_::Buffers_(mynest::iaf_cond_exp_sto_ip& n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

mynest::iaf_cond_exp_sto_ip::Buffers_::Buffers_(const Buffers_&, iaf_cond_exp_sto_ip& n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ---------------------------------------------------------------- 
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

mynest::iaf_cond_exp_sto_ip::iaf_cond_exp_sto_ip()
  : mynest::Archiving_Node_Sym(),
    P_(), 
    S_(P_),
    B_(*this)
{
  recordablesMap_.create();
}

mynest::iaf_cond_exp_sto_ip::iaf_cond_exp_sto_ip(const iaf_cond_exp_sto_ip& n)
  : mynest::Archiving_Node_Sym(n),
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{
}

mynest::iaf_cond_exp_sto_ip::~iaf_cond_exp_sto_ip()
{
  // GSL structs may not have been allocated, so we need to protect destruction
  if ( B_.s_ ) gsl_odeiv_step_free(B_.s_);
  if ( B_.c_ ) gsl_odeiv_control_free(B_.c_);
  if ( B_.e_ ) gsl_odeiv_evolve_free(B_.e_);
}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void mynest::iaf_cond_exp_sto_ip::init_state_(const Node& proto)
{
  const iaf_cond_exp_sto_ip& pr = downcast<iaf_cond_exp_sto_ip>(proto);
  S_ = pr.S_;
}

void mynest::iaf_cond_exp_sto_ip::init_buffers_()
{
  B_.spike_exc_.clear();          // includes resize
  B_.spike_inh_.clear();          // includes resize
  B_.currents_.clear();           // includes resize
  mynest::Archiving_Node_Sym::clear_history();

  B_.logger_.reset();

  B_.step_ = nest::Time::get_resolution().get_ms();
  B_.IntegrationStep_ = B_.step_;

  static const gsl_odeiv_step_type* T1 = gsl_odeiv_step_rkf45;
  
  if ( B_.s_ == 0 )
    B_.s_ = gsl_odeiv_step_alloc (T1, State_::STATE_VEC_SIZE);
  else 
    gsl_odeiv_step_reset(B_.s_);
    
  if ( B_.c_ == 0 )  
    B_.c_ = gsl_odeiv_control_y_new (1e-3, 0.0);
  else
    gsl_odeiv_control_init(B_.c_, 1e-3, 0.0, 1.0, 0.0);
    
  if ( B_.e_ == 0 )  
    B_.e_ = gsl_odeiv_evolve_alloc(State_::STATE_VEC_SIZE);
  else 
    gsl_odeiv_evolve_reset(B_.e_);
  
  B_.sys_.function  = iaf_cond_exp_sto_ip_dynamics;
  B_.sys_.jacobian  = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params    = reinterpret_cast<void*>(this);

  B_.I_stim_ = 0.0;
}

void mynest::iaf_cond_exp_sto_ip::calibrate()
{
  B_.logger_.init();  // ensures initialization in case mm connected after Simulate

  V_.RefractoryCounts_ = nest::Time(nest::Time::ms(P_.t_ref_abs)).get_steps();
  assert(V_.RefractoryCounts_ >= 0);  // since t_ref_ >= 0, this can only fail in error

  // Calculate inverse of the membrane capacitance
  V_.Time_step_in_s = nest::Time::get_resolution().get_ms()/1000.0;

  // Initialize random number generator
  V_.rng_ = net_->get_rng(get_thread());
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void mynest::iaf_cond_exp_sto_ip::update(nest::Time const & origin, const nest::long_t from, const nest::long_t to)
{
   
  assert(to >= 0 && (nest::delay) from < nest::Scheduler::get_min_delay());
  assert(from < to);

  for ( nest::long_t lag = from ; lag < to ; ++lag )
  {
    
    double t = 0.0;

    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t; this is of advantage
    // for a consistent and efficient integration across subsequent
    // simulation intervals
    while ( t < B_.step_ )
    {
      const int status = gsl_odeiv_evolve_apply(B_.e_, B_.c_, B_.s_, 
			   &B_.sys_,             // system of ODE
			   &t,                   // from t
			    B_.step_,            // to t <= step
			   &B_.IntegrationStep_, // integration step size
			    S_.y_); 	         // neuronal state

      if ( status != GSL_SUCCESS )
        throw nest::GSLSolverFailure(get_name(), status);
    }

    S_.y_[State_::G_EXC] += B_.spike_exc_.get_value(lag);
    S_.y_[State_::G_INH] += B_.spike_inh_.get_value(lag);
    
    // absolute refractory period
    if ( S_.r_ )
    {// neuron is absolute refractory
      --S_.r_; 
      S_.y_[State_::V_M] = P_.V_reset_; 
      S_.y_[State_::REFR] = 0.;
    }
    else
    {
      // Increment the time from the end of the absolute refractory period
      S_.time_rel_ref += B_.step_;

      // Calculate the refractoriness
      nest::double_t squared = S_.time_rel_ref * S_.time_rel_ref;
      S_.y_[State_::REFR] = squared / (P_.t_ref*P_.t_ref + squared);
      S_.y_[State_::GAIN] = S_.y_[State_::R_0] * log(1.+exp((S_.y_[State_::V_M] - S_.y_[State_::V_TH]) / S_.y_[State_::U_ALPHA]));
      S_.y_[State_::FIR_PROB] = S_.y_[State_::REFR] * S_.y_[State_::GAIN] * V_.Time_step_in_s;

      // Generate uniformly distributed random number in [0,1] range
      nest::double_t rand_num = V_.uni_dev_(V_.rng_);

      // neuron is decided to fire
      if ( rand_num <= S_.y_[State_::FIR_PROB] )
	  {
	      S_.r_              = V_.RefractoryCounts_;
	      S_.time_rel_ref	 = 0.0;
	      S_.y_[State_::V_M] = P_.V_reset_;

	      set_spiketime(nest::Time::step(origin.get_steps()+lag+1));
	  
	      nest::SpikeEvent se;
	      network()->send(*this, se, lag);
	  }
    }

    // Update intrinsic plasticity parameters
	nest::double_t old_gain = S_.y_[State_::R_0];
    nest::double_t old_thres = S_.y_[State_::V_TH];
    nest::double_t old_alpha = S_.y_[State_::U_ALPHA];
    //std::cout << "r0: " << old_gain << " Vth: " << old_thres << " Valpha: " << old_alpha << " FR: " << S_.y_[State_::GAIN] << " Vm: " << S_.y_[State_::V_M] << std::endl;
    S_.y_[State_::R_0] += P_.ip_rate * ( 1. - ( S_.y_[State_::GAIN] / P_.target_firing ) ) / old_gain;
    nest::double_t Aux = ( 1. + old_gain / P_.target_firing ) * (1. - exp(- S_.y_[State_::GAIN] / old_gain )) - 1.;
    S_.y_[State_::V_TH] += P_.ip_rate * Aux / old_alpha;
    S_.y_[State_::U_ALPHA] +=  P_.ip_rate * ( ( (S_.y_[State_::V_M] - old_thres) / old_alpha ) * Aux - 1. ) / old_alpha;

    // Check whether u_alpha is under 0
    //S_.y_[State_::U_ALPHA] = (S_.y_[State_::U_ALPHA]<=0)? std::numeric_limits<double_t>::min() : S_.y_[State_::U_ALPHA];

    // Check whether v_th is above 0
    //S_.y_[State_::V_TH] = (S_.y_[State_::V_TH]>=0)? 0.0 : S_.y_[State_::V_TH];


    // set new input current
    B_.I_stim_ = B_.currents_.get_value(lag);

    // log state data
    B_.logger_.record_data(origin.get_steps() + lag);

  }
}

void mynest::iaf_cond_exp_sto_ip::handle(nest::SpikeEvent & e)
{
  assert(e.get_delay() > 0);

  if(e.get_weight() > 0.0)
    B_.spike_exc_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			    e.get_weight() * e.get_multiplicity() );
  else
    B_.spike_inh_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			    -e.get_weight() * e.get_multiplicity() );  // ensure conductance is positive
}

void mynest::iaf_cond_exp_sto_ip::handle(nest::CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  const nest::double_t c=e.get_current();
  const nest::double_t w=e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
		      w *c);
}

void mynest::iaf_cond_exp_sto_ip::handle(nest::DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}

#endif //HAVE_GSL
